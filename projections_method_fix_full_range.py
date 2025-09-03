from itertools import combinations
from math import ceil, prod
import time
import tracemalloc
import numpy
from multiprocessing import Pool, cpu_count, Manager
from functools import partial
import psutil

#Вычисление ортогональных базисов
def calculate_orthogonal_bases(osn):
    full_diapazon = prod(osn)
    P = [full_diapazon // p for p in osn]
    beta = [P_i % p for P_i, p in zip(P, osn)]   
    m = []
    for i in range(len(osn)):
        try:
            m_val = pow(int(beta[i]), -1, int(osn[i]))
            m.append(m_val)
        except ValueError:
            return None
    
    result = [m_i * P_i for m_i, P_i in zip(m, P)]
    return result

#Вычисление скалярного произведения без переполнения
def safe_dot_product(a, b, mod):
    result = 0
    for a_i, b_i in zip(a, b):
        term = (a_i % mod) * (b_i % mod) % mod
        result = (result + term) % mod
    return result

#Решение системы уравнений
def solve_modular_system(M, b, mods):
    n = len(b)   
    if n == 1:
        try:
            inv = pow(int(M[0][0]), -1, int(mods[0]))
            solution = (b[0] * inv) % mods[0]
            return [solution]
        except Exception:
            return None
    
    M = [row[:] for row in M]
    b = b[:]
    for i in range(n):
        if M[i][i] == 0:
            for j in range(i + 1, n):
                if M[j][i] != 0:
                    M[i], M[j] = M[j], M[i]
                    b[i], b[j] = b[j], b[i]
                    break
        if M[i][i] == 0:
            return None
        try:
            inv = pow(int(M[i][i]), -1, int(mods[i]))
            for j in range(n):
                M[i][j] = (M[i][j] * inv) % mods[i]
            b[i] = (b[i] * inv) % mods[i]
        except Exception:
            return None
        for j in range(n):
            if j != i:
                factor = M[j][i]
                for k in range(n):
                    M[j][k] = (M[j][k] - factor * M[i][k]) % mods[j]
                b[j] = (b[j] - factor * b[i]) % mods[j]
    return b

#Вычисление методом проекций для одного числа
def check_single_number_with_projections(args):
    original_number, modular_number, osn_inform, osn_kontrol, osn, max_errors_mode, work_diapazon = args
    B = calculate_orthogonal_bases(osn)
    if B is None:
        return {
            'number': original_number,
            'modular': modular_number,
            'in_range': False,
            'A_value': 0,
            'corrections': None,
            'errors_fixed': 0,
            'corrected_A': 0
        }
    full_diapazon_value = prod(osn)
    A = safe_dot_product(modular_number, B, full_diapazon_value)
    in_range = A < work_diapazon
    result_info = {
        'number': original_number,
        'modular': modular_number,
        'in_range': in_range,
        'A_value': A,
        'corrections': None,
        'errors_fixed': 0,
        'corrected_A': 0
    }
    if not in_range:
        #Пытаемся исправить ошибку
        osn_kolvo = len(osn)
        k = len(osn_kontrol)
        if max_errors_mode == "гарантированно":
            max_errors_to_correct = k - 1
        else:
            max_errors_to_correct = k
        
        osn_int = [int(x) for x in osn]
        chisl_int = [int(x) for x in modular_number]
        for num_errors in range(1, max_errors_to_correct + 1):
            error_combinations = list(combinations(range(osn_kolvo), num_errors))
            for error_positions in error_combinations:
                proj_osn = [osn_int[i] for i in range(osn_kolvo) if i not in error_positions]
                proj_chisl = [chisl_int[i] for i in range(osn_kolvo) if i not in error_positions]
                B_proj = calculate_orthogonal_bases(proj_osn)
                if B_proj is None:
                    continue
                proj_full_diapazon = prod(proj_osn)
                A_proj = safe_dot_product(proj_chisl, B_proj, proj_full_diapazon)
                
                if A_proj < work_diapazon:
                    B_full = calculate_orthogonal_bases(osn_int)
                    if B_full is None:
                        continue
                    M = []
                    error_mods = [osn_int[pos] for pos in error_positions]
                    for i, pos_i in enumerate(error_positions):
                        row = []
                        for j, pos_j in enumerate(error_positions):
                            value = B_full[pos_i] % error_mods[j]
                            row.append(value)
                        M.append(row)
                    s = 0
                    for p in range(osn_kolvo):
                        if p not in error_positions:
                            term = (chisl_int[p] * B_full[p]) % full_diapazon_value
                            s = (s + term) % full_diapazon_value
                    b = []
                    for i, pos in enumerate(error_positions):
                        b_val = (A_proj - s) % osn_int[pos]
                        b.append(b_val)
                    try:
                        solutions = solve_modular_system(M, b, error_mods)
                        
                        if solutions is None:
                            continue
                        corrected = chisl_int[:]
                        error_positions_list = list(error_positions)
                        corrections_info = []
                        for i, pos in enumerate(error_positions_list):
                            old_value = corrected[pos]
                            corrected[pos] = solutions[i] % osn_int[pos]
                            corrections_info.append((pos, old_value, corrected[pos], osn_int[pos]))
                        A_corr = safe_dot_product(corrected, B_full, full_diapazon_value)
                        if A_corr < work_diapazon:
                            result_info['corrections'] = corrections_info
                            result_info['errors_fixed'] = num_errors
                            result_info['corrected_A'] = A_corr
                            return result_info
                    except Exception:
                        continue
    
    return result_info

#Генерация всех модулярных чисел
def generate_all_modular_numbers(osn):
    full_diapazon = prod(osn)
    numbers = []
    for i in range(int(full_diapazon)):
        a = [i % p for p in osn]
        numbers.append((i, a))
    return numbers

#Форматированный вывод информации о числе
def print_number_info(result_info, osn, work_diapazon):
    print(f"Проверка числа i = {result_info['number']}, a = {result_info['modular']}")    
    if result_info['in_range']:
        print(f"Число находится в пределах рабочего диапазона (0-{work_diapazon-1})")
        print(f"Вычисленное значение A = {result_info['A_value']}")
    else:
        print(f"Число находится ЗА ПРЕДЕЛАМИ рабочего диапазона (0-{work_diapazon-1})")
        print(f"Вычисленное значение A = {result_info['A_value']}")
        
        if result_info['corrections']:
            print(f"Приименен метод проекций для исправления ошибок")
            print(f"Обнаружено и исправлено {result_info['errors_fixed']} ошибок")
            print(f"Исправленное значение A = {result_info['corrected_A']}")
            
            print(f"\nИсправления:")
            for pos, old_value, new_value, modulus in result_info['corrections']:
                print(f"По основанию {modulus}: {old_value} → {new_value}")
            
            if result_info['corrected_A'] < work_diapazon:
                print(f"Исправление подтверждено: {result_info['corrected_A']} < {work_diapazon}")
            else:
                print(f"Исправление не прошло проверку: {result_info['corrected_A']} >= {work_diapazon}")
        else:
            print(f"Не удалось исправить ошибки методом проекций")

# Оценка размера файла
def estimate_file_size(full_diapazon, work_diapazon, osn):
    avg_line_length = len(f"Проверка числа i = 99999, a = {[999]*len(osn)}:\n") + 200
    total_size = full_diapazon * avg_line_length
    return total_size

def convert_bytes(size):
    for x in ['bytes', 'KB', 'MB', 'GB']:
        if size < 1024.0:
            return "%3.1f %s" % (size, x)
        size /= 1024.0

# Получение оптимального количества процессов
def get_optimal_process_count():
    try:
        #Используем psutil для получения информации о системе
        physical_cores = psutil.cpu_count(logical=False)
        logical_cores = psutil.cpu_count(logical=True)
        
        #Берем количество физических ядер + 50% для гипертрединга
        optimal_processes = min(physical_cores * 1.5, logical_cores)
        return max(1, int(optimal_processes))

    except:
        #Если psutil недоступен, используем безопасное значение
        return min(cpu_count(), 8)  # Ограничиваем 8 процессами по умолчанию

#Основная функция проверки с многопроцессорностью
def check_all_numbers_parallel(osn_inform, osn_kontrol, max_errors_mode="гарантированно", verbose=False):
    osn = osn_inform + osn_kontrol
    work_diapazon = prod(osn_inform)
    full_diapazon = prod(osn)
    
    print(f"Информация о системе:")
    print(f"Информационные основания: {osn_inform}")
    print(f"Контрольные основания: {osn_kontrol}")
    print(f"Все основания: {osn}")
    print(f"Рабочий диапазон: 0-{work_diapazon-1}")
    print(f"Полный диапазон: 0-{full_diapazon-1}")
    print(f"Количество чисел: {full_diapazon}")
    print(f"Режим исправления: {max_errors_mode}")
    #Получаем оптимальное количество процессов
    num_processes = get_optimal_process_count()
    print(f"Оптимальное количество процессов: {num_processes}")
    #Генерируем все числа
    numbers = generate_all_modular_numbers(osn)
    #Подготавливаем аргументы для параллельной обработки
    task_args = [(num, mod, osn_inform, osn_kontrol, osn, max_errors_mode, work_diapazon) 
                 for num, mod in numbers]
    
    print(f"Запуск параллельной обработки с {num_processes} процессами")
    
    stats = {
        'total': 0,
        'correct': 0,
        'errors_fixed': {i: 0 for i in range(1, len(osn_kontrol) + 1)},
        'errors_not_fixed': 0,
        'detailed_results': []
    }
    
    #Настройка вывода на экран статуса выполнения (каждые 5%)
    progress_step = max(1, full_diapazon // 20)  #5% от общего количества
    next_progress_threshold = progress_step
    with Pool(processes=num_processes) as pool:
        results = pool.imap_unordered(check_single_number_with_projections, task_args)
        for i, result_info in enumerate(results, 1):
            stats['total'] = i
            stats['detailed_results'].append(result_info)
            if result_info['in_range']:
                stats['correct'] += 1
            elif result_info['corrections'] and result_info['corrected_A'] < work_diapazon:
                stats['errors_fixed'][result_info['errors_fixed']] += 1
            else:
                stats['errors_not_fixed'] += 1
            if verbose:
                print_number_info(result_info, osn, work_diapazon)
            #Вывод прогресса каждые 5%
            if i >= next_progress_threshold:
                percentage = (i / full_diapazon) * 100
                print(f"Обработано: {i}/{full_diapazon} ({percentage:.1f}%)")
                next_progress_threshold += progress_step
            #Также выводим каждые 10000 чисел для больших систем
            if full_diapazon > 100000 and i % 10000 == 0:
                percentage = (i / full_diapazon) * 100
                print(f"Промежуточный прогресс: {i}/{full_diapazon} ({percentage:.1f}%)")
    
    return stats

#Запись результата в файл
def save_results_to_file(stats, osn_inform, osn_kontrol, max_errors_mode, filename):
    osn = osn_inform + osn_kontrol
    work_diapazon = prod(osn_inform)
    
    with open(filename, 'w', encoding='utf-8') as file:
        file.write("Полные результаты проверки модулярных чисел\n")
        file.write("Информация о системе:\n")
        file.write(f"Информационные основания: {osn_inform}\n")
        file.write(f"Контрольные основания: {osn_kontrol}\n")
        file.write(f"Все основания: {osn}\n")
        file.write(f"Рабочий диапазон: 0-{work_diapazon-1}\n")
        file.write(f"Полный диапазон: 0-{prod(osn)-1}\n")
        file.write(f"Режим исправления: {max_errors_mode}\n\n")
        file.write("Итоговая статистика проверки\n")
        file.write(f"Всего чисел проверено: {stats['total']}\n")
        file.write(f"Корректных чисел: {stats['correct']}\n")
        file.write(f"Чисел с ошибками, которые удалось исправить: {sum(stats['errors_fixed'].values())}\n")
        file.write(f"Чисел с ошибками, которые не удалось исправить: {stats['errors_not_fixed']}\n\n")
        file.write("Исправленные ошибки по кратности:\n")
        for num_errors, count in stats['errors_fixed'].items():
            if count > 0:
                file.write(f"Ошибок кратности {num_errors}: {count}\n")
        
        file.write("Подробные результаты по каждому числу\n")
        
        for result in stats['detailed_results']:
            file.write(f"Число i = {result['number']}, a = {result['modular']}\n")
            if result['in_range']:
                file.write(f"В пределах рабочего диапазона, A = {result['A_value']}\n")
            else:
                file.write(f"За пределами рабочего диапазона, A = {result['A_value']}\n")
                if result['corrections']:
                    file.write(f"Исправлено {result['errors_fixed']} ошибок, исправленное A = {result['corrected_A']}\n")
                    for pos, old_value, new_value, modulus in result['corrections']:
                        file.write(f"По основанию {modulus}: {old_value} → {new_value}\n")
                else:
                    file.write("Не удалось исправить ошибки\n")
                    
        file.write("Ресурсы программы\n")
        current, peak = tracemalloc.get_traced_memory()
        execution_time = time.time() - start_time
        file.write(f"Пиковое использование памяти: {peak / 1024:.2f} KB\n")
        file.write(f"Время выполнения: {execution_time:.6f} секунд\n")

def main():
    global start_time
    start_time = time.time()
    tracemalloc.start()
    #Выбор режима работы
    print("\nВыберите режим исправления ошибок:")
    print("1 - Гарантированное исправление (до k-1 ошибок)")
    print("2 - Попытка исправления (до k ошибок, не гарантировано)")
    mode_choice = input("Ваш выбор (1 или 2): ").strip()
    
    if mode_choice == "1":
        max_errors_mode = "гарантированно"
    else:
        max_errors_mode = "попытка"
    
    #Ввод данных
    osn_inform_kolvo = int(input('\nВведите количество информационных оснований: '))
    osn_kontrol_kolvo = int(input('Введите количество контрольных оснований: '))
    
    osn_inform = []
    osn_kontrol = []
    
    for i in range(osn_inform_kolvo):
        osn_inform.append(int(input(f'Введите {i+1} информационное основание: ')))
    
    for i in range(osn_kontrol_kolvo):
        osn_kontrol.append(int(input(f'Введите {i+1} контрольное основание: ')))
    
    osn = osn_inform + osn_kontrol
    full_diapazon = prod(osn)
    
    #Проверяем все числа с использованием многопроцессорности
    print(f"\nЗапуск проверки {full_diapazon} чисел...")
    stats = check_all_numbers_parallel(osn_inform, osn_kontrol, max_errors_mode, verbose=False)
    
    #Вывод статистики
    print("Итоговая статистика проверки")
    print(f"Всего чисел проверено: {stats['total']}")
    print(f"Корректных чисел: {stats['correct']}")
    print(f"Чисел с ошибками, которые удалось исправить: {sum(stats['errors_fixed'].values())}")
    print(f"Чисел с ошибками, которые не удалось исправить: {stats['errors_not_fixed']}")
    print("\nИсправленные ошибки по кратности:")
    for num_errors, count in stats['errors_fixed'].items():
        if count > 0:
            print(f"  Ошибок кратности {num_errors}: {count}")
    #Замер ресурсов
    current, peak = tracemalloc.get_traced_memory()
    execution_time = time.time() - start_time
    print(f"\nРесурсы программы:")
    print(f"Пиковое использование памяти: {peak / 1024:.2f} KB")
    print(f"Время выполнения: {execution_time:.6f} секунд")
    
    #Предложение сохранить результаты в файл
    estimated_size = estimate_file_size(full_diapazon, prod(osn_inform), osn)
    print(f"\nПриблизительный размер файла с полными результатами: {convert_bytes(estimated_size)}")
    save_choice = input("\nЗаписать полные результаты в файл? [yes/no]: ").strip().lower()
    if save_choice == 'yes':
        filename = input("Введите имя файла (по умолчанию: results_parallel.txt): ").strip()
        if not filename:
            filename = "results_parallel.txt"
        save_results_to_file(stats, osn_inform, osn_kontrol, max_errors_mode, filename)
        print(f"Результаты сохранены в файл: {filename}")
    
    #Предложение вывести результаты на экран
    display_choice = input("\nВывести подробные результаты на экран? [yes/no]: ").strip().lower()
    if display_choice == 'yes':
        print("Подробные результаты проверки")
        
        for result in stats['detailed_results']:
            print_number_info(result, osn, prod(osn_inform))

main()
