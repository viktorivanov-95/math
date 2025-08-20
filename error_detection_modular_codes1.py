import numpy as np
import sys
from collections import defaultdict
from itertools import combinations, product

def calculate_distances(code, codes_to_compare):
    distances = []
    for other_code in codes_to_compare:
        if np.array_equal(code, other_code):
            continue  
        dist = np.sum(code != other_code)
        distances.append(dist)
    
    if not distances:
        return 0, 0, []
    
    return min(distances), max(distances), distances

def analyze_error_correction(code, inform_codes, osn):
    n = len(osn)
    error_stats = {i: {'total': 0, 'correctable': 0} for i in range(1, n+1)}
    
    for error_count in range(1, n+1):
        for error_positions in combinations(range(n), error_count):
            error_values = []
            for pos in error_positions:
                error_values.append([val for val in range(osn[pos]) if val != code[pos]])
            
            for error_combination in product(*error_values):
                erroneous_code = code.copy()
                for i, pos in enumerate(error_positions):
                    erroneous_code[pos] = error_combination[i]
                
                is_in_inform_range = any(np.array_equal(erroneous_code, c) for c in inform_codes)
                
                error_stats[error_count]['total'] += 1
                if not is_in_inform_range:
                    error_stats[error_count]['correctable'] += 1
    
    result = {}
    for error_count in range(1, n+1):
        if error_stats[error_count]['total'] > 0:
            percentage = (error_stats[error_count]['correctable'] / error_stats[error_count]['total']) * 100
            result[error_count] = {
                'percentage': percentage,
                'correctable': error_stats[error_count]['correctable'],
                'total': error_stats[error_count]['total']
            }
        else:
            result[error_count] = {
                'percentage': 0,
                'correctable': 0,
                'total': 0
            }
    
    return result

def estimate_file_size(full_diapazon, work_diapazon, osn):
    avg_inform_line = len(f"i = 99999, a = {np.array([999]*len(osn))}, min_dist = 99, max_dist = 99, 1 ошибок - 100.00% (9999/9999), 2 ошибок - 100.00% (9999/9999), 3 ошибок - 99.99% (9999/9999), 4 ошибок - 99.99% (9999/9999), 5 ошибок - 99.99% (9999/9999)\n")
    header_lines = 15
    avg_header_line = 50
    
    total_size = (work_diapazon * avg_inform_line + header_lines * avg_header_line)
    
    return total_size

def convert_bytes(size):
    for x in ['bytes', 'KB', 'MB', 'GB']:
        if size < 1024.0:
            return "%3.1f %s" % (size, x)
        size /= 1024.0

def main():
    osn_inform_kolvo = int(input('Введите количество информационных оснований: '))
    osn_kontrol_kolvo = int(input('Введите количество контрольных оснований: '))
    
    osn_inform = []
    osn_kontrol = []
    
    for i in range(osn_inform_kolvo):
        osn_inform.append(int(input(f'Введите {i+1} информационное основание: ')))
    for i in range(osn_kontrol_kolvo):
        osn_kontrol.append(int(input(f'Введите {i+1} контрольное основание: ')))
    
    osn = osn_inform + osn_kontrol
    
    work_diapazon = 1
    for base in osn_inform:
        work_diapazon *= base
    
    full_diapazon = work_diapazon
    for base in osn_kontrol:
        full_diapazon *= base
    
    print(f'\nКоличество чисел в разрешенном диапазоне (0-{work_diapazon-1}): {work_diapazon}')
    print(f'Количество чисел в запрещенном диапазоне ({work_diapazon}-{full_diapazon-1}): {full_diapazon - work_diapazon}')
    print(f'Количество всех чисел (0-{full_diapazon-1}): {full_diapazon}')
    
    estimated_size = estimate_file_size(full_diapazon, work_diapazon, osn)
    print(f"\nПриблизительный размер файла с результатами: {convert_bytes(estimated_size)}")
    
    filename = 'modular_numbers_output.txt'
    save_to_file = input(f'Записать все числа с расстояниями в файл {filename}? [yes/no]: ').strip().lower()
    
    file_created = False
    inform_min_dists = []
    inform_max_dists = []
    
    # Для статистики по ошибкам
    error_stats_sum = {i: {'correctable': 0, 'total': 0} for i in range(1, len(osn)+1)}
    
    print("\nВычисление информационных кодов...")
    inform_codes = []
    for i in range(work_diapazon):
        code = np.array([i % base for base in osn], dtype=np.int64)
        inform_codes.append(code)
        if i % 100000 == 0 and i > 0:
            print(f"Обработано {i}/{work_diapazon} информационных кодов")
    
    if save_to_file == 'yes':
        with open(filename, 'w', encoding='utf-8') as file:
            file_created = True
            
            file.write(f'Информационные основания: {osn_inform}\n')
            file.write(f'Контрольные основания: {osn_kontrol}\n')
            file.write(f'Все основания: {osn}\n')
            file.write(f'Рабочий диапазон: 0-{work_diapazon-1}\n')
            file.write(f'Полный диапазон: 0-{full_diapazon-1}\n\n')
            file.write(f'Количество чисел в разрешенном диапазоне: {work_diapazon}\n')
            file.write(f'Количество чисел в запрещенном диапазоне: {full_diapazon - work_diapazon}\n')
            file.write(f'Количество всех чисел: {full_diapazon}\n\n')
            
            file.write('===== Статистика по расстояниям и исправлению ошибок =====\n')
            
            overall_inform_min = float('inf')
            overall_inform_max = 0
            
            file.write('\n===== Разрешенный диапазон (информационный) с расстояниями и анализом ошибок =====\n')
            print("\nВычисление расстояний и анализ ошибок для информационного диапазона...")
            
            for i in range(work_diapazon):
                code = inform_codes[i]
                min_dist, max_dist, dist_array = calculate_distances(code, inform_codes)
                error_analysis = analyze_error_correction(code, inform_codes, osn)
                
                inform_min_dists.append(min_dist)
                inform_max_dists.append(max_dist)
                
                # Обновляем статистику по ошибкам
                for error_count in error_analysis:
                    error_stats_sum[error_count]['correctable'] += error_analysis[error_count]['correctable']
                    error_stats_sum[error_count]['total'] += error_analysis[error_count]['total']
                
                # Формируем строку с анализом ошибок
                error_parts = []
                for error_count in sorted(error_analysis.keys()):
                    stats = error_analysis[error_count]
                    error_parts.append(f"{error_count} ошибок - {stats['percentage']:.2f}% ({stats['correctable']}/{stats['total']})")
                
                error_str = ", ".join(error_parts)
                
                file.write(f'i = {i}, a = {code}, min_dist = {min_dist}, max_dist = {max_dist}, {error_str}\n')
                
                if min_dist < overall_inform_min:
                    overall_inform_min = min_dist
                if max_dist > overall_inform_max:
                    overall_inform_max = max_dist
                
                if i % 100 == 0 and i > 0:
                    print(f"Обработано {i}/{work_diapazon} информационных чисел")
            
            # Запись общей статистики
            file.write('\n===== Общая статистика =====\n')
            file.write(f'Минимальное расстояние в разрешенном диапазоне: {overall_inform_min}\n')
            file.write(f'Максимальное расстояние в разрешенном диапазоне: {overall_inform_max}\n')
            
            # Статистика по исправлению ошибок
            file.write('\nСредняя эффективность обнаружения ошибок:\n')
            for error_count in sorted(error_stats_sum.keys()):
                if error_stats_sum[error_count]['total'] > 0:
                    avg_percentage = (error_stats_sum[error_count]['correctable'] / error_stats_sum[error_count]['total']) * 100
                    file.write(f'{error_count} ошибок: {avg_percentage:.2f}% '
                              f'({error_stats_sum[error_count]["correctable"]}/{error_stats_sum[error_count]["total"]})\n')
                else:
                    file.write(f'{error_count} ошибок: 0.00% (0/0)\n')
        
        print(f'Файл успешно создан: {filename}')
    else:
        print("\nВычисление расстояний для информационного диапазона...")
        for i in range(work_diapazon):
            code = inform_codes[i]
            min_dist, max_dist, _ = calculate_distances(code, inform_codes)
            inform_min_dists.append(min_dist)
            inform_max_dists.append(max_dist)
            
            # Все равно собираем статистику по ошибкам
            error_analysis = analyze_error_correction(code, inform_codes, osn)
            for error_count in error_analysis:
                error_stats_sum[error_count]['correctable'] += error_analysis[error_count]['correctable']
                error_stats_sum[error_count]['total'] += error_analysis[error_count]['total']
            
            if i % 100000 == 0 and i > 0:
                print(f"Обработано {i}/{work_diapazon} информационных чисел")
        
        print('Запись в файл пропущена.')
    
    # Вывод статистики на экран
    print('\n===== Статистика по расстояниям =====')
    print(f'Минимальное расстояние в разрешенном диапазоне: {min(inform_min_dists)}')
    print(f'Максимальное расстояние в разрешенном диапазоне: {max(inform_max_dists)}')
    
    # Вывод средней статистики по ошибкам
    print('\n===== Средняя эффективность обнаружения ошибок =====')
    for error_count in sorted(error_stats_sum.keys()):
        if error_stats_sum[error_count]['total'] > 0:
            avg_percentage = (error_stats_sum[error_count]['correctable'] / error_stats_sum[error_count]['total']) * 100
            print(f'{error_count} ошибок: {avg_percentage:.2f}% '
                  f'({error_stats_sum[error_count]["correctable"]}/{error_stats_sum[error_count]["total"]})')
        else:
            print(f'{error_count} ошибок: 0.00% (0/0)')
    
    # Анализ ошибок для нескольких примеров
    print('\n===== Анализ исправления ошибок (первые 3 числа) =====')
    for i in range(min(3, work_diapazon)):
        code = inform_codes[i]
        error_analysis = analyze_error_correction(code, inform_codes, osn)
        error_parts = []
        for error_count in sorted(error_analysis.keys()):
            stats = error_analysis[error_count]
            error_parts.append(f"{error_count} ошибок - {stats['percentage']:.2f}% ({stats['correctable']}/{stats['total']})")
        error_str = ", ".join(error_parts)
        print(f'i = {i}, a = {code}, {error_str}')

if __name__ == "__main__":
    main()
