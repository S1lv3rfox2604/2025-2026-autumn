#include <cmath>
#include <algorithm>
#include <iostream>
#include <queue>
#include <pthread.h>
#include "solver.h"

// Структуры для пула потоков
typedef struct {
    int type; // 0 - поиск максимума, 1 - исключение
    int k;    // текущий шаг
    int start_row;
    int end_row;
} task_t;

typedef struct {
    int n;
    double* A;
    double* b;
    int* col_perm;
    pthread_mutex_t queue_mutex;
    pthread_cond_t queue_cond;
    pthread_cond_t tasks_done_cond;
    std::queue<task_t>* task_queue;
    bool shutdown;
    int tasks_completed;
    int total_tasks;
    
    // Для поиска максимума
    pthread_mutex_t max_mutex;
    double global_max_val;
    int global_i_max;
    int global_j_max;
    
    bool matrix_singular;
} thread_pool_t;

static thread_pool_t pool;
static pthread_t* threads = nullptr;

// Функция рабочего потока
void* worker_thread(void* arg) {
    int thread_id = *(int*)arg;
    
    while (true) {
        pthread_mutex_lock(&pool.queue_mutex);
        
        // Ожидание задачи
        while (pool.task_queue->empty() && !pool.shutdown) {
            pthread_cond_wait(&pool.queue_cond, &pool.queue_mutex);
        }
        
        // Завершение работы
        if (pool.shutdown && pool.task_queue->empty()) {
            pthread_mutex_unlock(&pool.queue_mutex);
            break;
        }
        
        // Получение задачи
        task_t task = pool.task_queue->front();
        pool.task_queue->pop();
        pthread_mutex_unlock(&pool.queue_mutex);
        
        // Выполнение задачи
        if (task.type == 0) { // Поиск максимума
            double local_max_val = 0.0;
            int local_i_max = task.k;
            int local_j_max = task.k;
            
            for (int i = task.start_row; i < task.end_row; ++i) {
                if (i < pool.n) {
                    double* row_i = pool.A + i * pool.n;
                    for (int j = task.k; j < pool.n; ++j) {
                        double val = std::abs(row_i[pool.col_perm[j]]);
                        if (val > local_max_val) {
                            local_max_val = val;
                            local_i_max = i;
                            local_j_max = j;
                        }
                    }
                }
            }
            
            // Обновление глобального максимума
            pthread_mutex_lock(&pool.max_mutex);
            if (local_max_val > pool.global_max_val) {
                pool.global_max_val = local_max_val;
                pool.global_i_max = local_i_max;
                pool.global_j_max = local_j_max;
            }
            pthread_mutex_unlock(&pool.max_mutex);
            
        } else if (task.type == 1) { // Зануление
            double* row_k = pool.A + task.k * pool.n;
            double pivot = row_k[pool.col_perm[task.k]];
            
            for (int i = task.start_row; i < task.end_row; ++i) {
                if (i < pool.n && i > task.k) { // Только строки ниже диагонали
                    double* row_i = pool.A + i * pool.n;
                    double factor = row_i[pool.col_perm[task.k]] / pivot;
                    
                    for (int j = task.k; j < pool.n; ++j) {
                        row_i[pool.col_perm[j]] -= factor * row_k[pool.col_perm[j]];
                    }
                    pool.b[i] -= factor * pool.b[task.k];
                }
            }
        }
        
        // Уведомление о завершении задачи
        pthread_mutex_lock(&pool.queue_mutex);
        pool.tasks_completed++;
        if (pool.tasks_completed == pool.total_tasks) {
            pthread_cond_signal(&pool.tasks_done_cond);
        }
        pthread_mutex_unlock(&pool.queue_mutex);
    }
    
    pthread_exit(NULL);
}

// Функция добавления задач и ожидания их завершения
void execute_tasks(int type, int k, int num_tasks) {
    // Сброс счетчиков
    pthread_mutex_lock(&pool.queue_mutex);
    pool.tasks_completed = 0;
    pool.total_tasks = num_tasks;
    pthread_mutex_unlock(&pool.queue_mutex);
    
    // Добавление задач
    int rows_per_task = (pool.n - ((type == 1) ? k + 1 : k)) / num_tasks;
    int extra_rows = (pool.n - ((type == 1) ? k + 1 : k)) % num_tasks;
    
    int current_row = (type == 1) ? k + 1 : k;
    
    for (int i = 0; i < num_tasks; ++i) {
        task_t task;
        task.type = type;
        task.k = k;
        task.start_row = current_row;
        
        int rows_for_this_task = rows_per_task + (i < extra_rows ? 1 : 0);
        task.end_row = current_row + rows_for_this_task;
        
        // Корректировка для последней задачи
        if (i == num_tasks - 1) {
            task.end_row = (type == 1) ? pool.n : pool.n;
        }
        
        current_row = task.end_row;
        
        pthread_mutex_lock(&pool.queue_mutex);
        pool.task_queue->push(task);
        pthread_cond_signal(&pool.queue_cond);
        pthread_mutex_unlock(&pool.queue_mutex);
    }
    
    // Ожидание завершения всех задач
    pthread_mutex_lock(&pool.queue_mutex);
    while (pool.tasks_completed < pool.total_tasks) {
        pthread_cond_wait(&pool.tasks_done_cond, &pool.queue_mutex);
    }
    pthread_mutex_unlock(&pool.queue_mutex);
}

// Параллельная версия с пулом потоков
bool gauss_solve_parallel(int n, double* A, double* b, double* x, int num_threads) {
    int* col_perm = new int[n];
    for (int j = 0; j < n; ++j) col_perm[j] = j;
    
    // Инициализация пула
    pool.n = n;
    pool.A = A;
    pool.b = b;
    pool.col_perm = col_perm;
    pool.shutdown = false;
    pool.matrix_singular = false;
    pool.task_queue = new std::queue<task_t>();
    pool.tasks_completed = 0;
    pool.total_tasks = 0;
    
    pthread_mutex_init(&pool.queue_mutex, NULL);
    pthread_mutex_init(&pool.max_mutex, NULL);
    pthread_cond_init(&pool.queue_cond, NULL);
    pthread_cond_init(&pool.tasks_done_cond, NULL);
    
    // Создание рабочих потоков
    threads = new pthread_t[num_threads];
    int* thread_ids = new int[num_threads];
    
    for (int i = 0; i < num_threads; ++i) {
        thread_ids[i] = i;
        pthread_create(&threads[i], NULL, worker_thread, &thread_ids[i]);
    }
    
    bool success = true;
    
    // Прямой ход метода Гаусса
    for (int k = 0; k < n; ++k) {
        // Поиск главного элемента
        pool.global_max_val = 0.0;
        pool.global_i_max = k;
        pool.global_j_max = k;
        
        execute_tasks(0, k, num_threads);
        
        // Проверка на вырожденность
        if (pool.global_max_val == 0.0) {
            success = false;
            break;
        }
        
        // Перестановка строк
        if (pool.global_i_max != k) {
            double* row_k = A + k * n;
            double* row_im = A + pool.global_i_max * n;
            for (int j = 0; j < n; ++j) {
                std::swap(row_k[j], row_im[j]);
            }
            std::swap(b[k], b[pool.global_i_max]);
        }
        
        // Перестановка столбцов
        if (pool.global_j_max != k) {
            std::swap(col_perm[k], col_perm[pool.global_j_max]);
        }
        
        // Исключение переменной
        if (k < n - 1) {
            execute_tasks(1, k, num_threads);
        }
    }
    
    // Завершение работы пула
    pthread_mutex_lock(&pool.queue_mutex);
    pool.shutdown = true;
    pthread_cond_broadcast(&pool.queue_cond);
    pthread_mutex_unlock(&pool.queue_mutex);
    
    // Ожидание завершения потоков
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], NULL);
    }
    
    // Обратный ход
    if (success) {
        for (int i = n-1; i >= 0; --i) {
            double* row_i = A + i * n;
            double sum = b[i];
            for (int j = i+1; j < n; ++j) {
                sum -= row_i[col_perm[j]] * x[j];
            }
            if (std::abs(row_i[col_perm[i]]) < 1e-15) {
                success = false;
                break;
            }
            x[i] = sum / row_i[col_perm[i]];
        }
    }
    
    // Освобождение ресурсов
    delete[] threads;
    delete[] thread_ids;
    delete pool.task_queue;
    
    pthread_mutex_destroy(&pool.queue_mutex);
    pthread_mutex_destroy(&pool.max_mutex);
    pthread_cond_destroy(&pool.queue_cond);
    pthread_cond_destroy(&pool.tasks_done_cond);
    
    delete[] col_perm;
    return success;
}