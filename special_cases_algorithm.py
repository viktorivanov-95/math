import numpy as np
import math
from multiprocessing import Pool, cpu_count
import psutil
import time
from collections import defaultdict

#Находим минимальное расстояние и соответствующее число из информационного диапазона
def find_min_distance_info(code, inform_codes, inform_numbers):
    min_dist = float('inf')
    min_dist_info_number = None
    min_dist_info_code = None
    for other_code, other_number in zip(inform_codes, inform_numbers):
        dist = np.sum(code != other_code)
        if dist < min_dist:
            min_dist = dist
            min_dist_info_number = other_number
            min_dist_info_code = other_code
    return min_dist, min_dist_info_number, min_dist_info_code

#Находим все информационные коды на минимальном расстоянии
def find_all_min_distance_codes(code, inform_codes, inform_numbers, min_dist):
    min_dist_codes = []
    min_dist_numbers = []
    for other_code, other_number in zip(inform_codes, inform_numbers):
        dist = np.sum(code != other_code)
        if dist == min_dist:
            min_dist_codes.append(other_code)
            min_dist_numbers.append(other_number)
    return min_dist_codes, min_dist_numbers

#Генерируем модулярное представление числа A
def generate_modular_number(i, osn):
    return np.array([i % p for p in osn])

#Обрабатываем чанк чисел с детальным анализом расстояний
def process_chunk_detailed(args):
    chunk_start, chunk_end, osn, inform_codes, inform_numbers, total_mods = args
    results = []
    for i in range(chunk_start, chunk_end):
        code = generate_modular_number(i, osn)
        #Находим минимальное расстояние и соответствующее информационное число
        min_dist, min_info_number, min_info_code = find_min_distance_info(code, inform_codes, inform_numbers)
        #Находим все информационные коды на минимальном расстоянии
        all_min_codes, all_min_numbers = find_all_min_distance_codes(code, inform_codes, inform_numbers, min_dist)
        #Вычисляем все расстояния до информационных кодов
        all_distances = []
        for other_code in inform_codes:
            dist = np.sum(code != other_code)
            all_distances.append(dist)
        #Считаем распределение расстояний ТОЛЬКО для этого числа
        dist_distribution = {}
        for dist in range(1, total_mods + 1):  # Все расстояния от 1 до n
            count = all_distances.count(dist)
            dist_distribution[dist] = count
        results.append((i, code, min_dist, min_info_number, min_info_code, 
                       all_min_codes, all_min_numbers, dist_distribution, all_distances))
    return results
  
#Определяем оптимальное количество процессов
def get_optimal_process_count():
    try:
        physical_cores = psutil.cpu_count(logical=False)
        logical_cores = psutil.cpu_count(logical=True)
        optimal_processes = min(physical_cores * 1.5, logical_cores)
        return max(1, int(optimal_processes))
    except:
        return min(cpu_count(), 8)

#Сохранение особых случаев в файл
def save_special_cases(filename, special_cases, osn_inform):
    with open(filename, 'w', encoding='utf-8') as f:
        f.write("Особые случаи: числа с расстояниями больше их min_dist но меньше или равно заданному пользователем max_min_dist\n")        
        for case in special_cases:
            (number, code, min_dist, min_info_number, min_info_code, 
             special_dist, special_info_number, special_info_code) = case
            
            f.write(f"Число из запрещенного диапазона:\n")
            f.write(f"  i = {number}, a = {code}\n")
            f.write(f"  min_dist = {min_dist}\n")
            f.write(f"  Ближайшее информационное число: i = {min_info_number}, a = {min_info_code}\n")
            f.write(f"\nОсобое расстояние:\n")
            f.write(f"  Расстояние {special_dist} до информационного числа: i = {special_info_number}, a = {special_info_code}\n")
            f.write("-" * 80 + "\n\n")

#Сохраняем случаи с несколькими кодами на min_dist
def save_multiple_min_cases(filename, multiple_min_cases, osn_inform):
    with open(filename, 'w', encoding='utf-8') as f:
        f.write("Особые случаи: числа с несколькими информационными кодами на минимальном расстоянии\n")        
        for case in multiple_min_cases:
            number, code, min_dist, all_min_numbers, all_min_codes = case
            f.write(f"Число из запрещенного диапазона:\n")
            f.write(f"  i = {number}, a = {code}\n")
            f.write(f"  min_dist = {min_dist}\n")
            f.write(f"  Количество информационных кодов на min_dist: {len(all_min_numbers)}\n")
            f.write(f"  Все информационные числа на min_dist: {all_min_numbers}\n")
            for i, (info_num, info_code) in enumerate(zip(all_min_numbers, all_min_codes)):
                f.write(f"  Вариант {i+1}: i = {info_num}, a = {info_code}\n")
            f.write("-" * 80 + "\n\n")

#Получает группы min_dist
def get_user_min_dist_groups(total_mods):
    print("НАСТРОЙКА АНАЛИЗИРУЕМЫХ ГРУПП MIN_DIST")    
    print(f"\nВсего возможных расстояний: от 1 до {total_mods}")
    print("Вы можете задать максимальное значение min_dist для анализа.")
    print("Будут показаны все группы от 1 до указанного значения.")
    while True:
        try:
            max_min_dist = int(input(f"\nВведите максимальное значение min_dist для анализа (1-{total_mods}): "))
            if 1 <= max_min_dist <= total_mods:
                return list(range(1, max_min_dist + 1)), max_min_dist
            else:
                print(f"Ошибка: введите число от 1 до {total_mods}")
        except ValueError:
            print("Ошибка: введите целое число")

#Основная функция анализа модулярного кода с пользовательскими группами min_dist
def analyze_modular_code_final():
    #Ввод данных
    print("Анализ возможностей модулярного кода по обнаружению и исправлению ошибок")    
    osn_inform_kolvo = int(input('Введите количество информационных оснований: '))
    osn_kontrol_kolvo = int(input('Введите количество контрольных оснований: '))
    osn_inform = []
    osn_kontrol = []
    for i in range(osn_inform_kolvo):
        osn_inform.append(int(input(f'Введите {i+1} информационное основание: ')))
    for i in range(osn_kontrol_kolvo):
        osn_kontrol.append(int(input(f'Введите {i+1} контрольное основание: ')))
    #Объединение оснований
    osn = osn_inform + osn_kontrol
    total_mods = len(osn)
    #Получение интересующих групп min_dist от пользователя
    min_dist_groups, user_max_min_dist = get_user_min_dist_groups(total_mods)
    #Вычисление диапазонов
    work_diapazon = math.prod(osn_inform)
    full_diapazon = math.prod(osn)
    print(f'\nКоличество чисел в разрешенном диапазоне (0-{work_diapazon-1}): {work_diapazon:,}')
    print(f'Количество чисел в запрещенном диапазоне ({work_diapazon}-{full_diapazon-1}): {full_diapazon - work_diapazon:,}')
    print(f'Количество всех чисел (0-{full_diapazon-1}): {full_diapazon:,}')
    #Теоретические возможности кода
    print(f"\nТеоретические возможности кода (по теории модулярного кодирования):")
    print(f"Количество контрольных оснований: {osn_kontrol_kolvo}")
    print(f"Гарантированно обнаруживаемые ошибки: кратность 1 - {osn_kontrol_kolvo}")
    if osn_kontrol_kolvo >= 2:
        print(f"Гарантированно исправляемые ошибки: кратность 1 - {osn_kontrol_kolvo - 1}")
    #Генерация информационных кодов
    print("\nГенерация информационных кодов...")
    inform_codes = []
    inform_numbers = []  #Сохраняем также исходные числа
    for i in range(work_diapazon):
        inform_codes.append(generate_modular_number(i, osn))
        inform_numbers.append(i)
    #Подготовка к многопроцессорной обработке
    num_processes = get_optimal_process_count()
    print(f"Используется процессов: {num_processes}")
    #Разделение запрещенного диапазона на чанки
    chunk_size = (full_diapazon - work_diapazon) // num_processes
    chunks = []
    
    for i in range(num_processes):
        chunk_start = work_diapazon + i * chunk_size
        chunk_end = chunk_start + chunk_size
        if i == num_processes - 1:  #Последний чанк включает остаток
            chunk_end = full_diapazon
        chunks.append((chunk_start, chunk_end, osn, inform_codes, inform_numbers, total_mods))
    
    #Многопроцессорная обработка
    print("Анализ минимальных расстояний для запрещенного диапазона...")
    start_time = time.time()
    with Pool(processes=num_processes) as pool:
        results = pool.map(process_chunk_detailed, chunks)
    
    #Сбор и анализ результатов
    min_dist_counts = {}
    group_distributions = defaultdict(lambda: defaultdict(int))
    total_distances_per_group = defaultdict(int)
    special_cases = []  #Для хранения особых случаев
    multiple_min_cases = []  #Для хранения случаев с несколькими кодами на min_dist
    #Статистика по множественным минимальным расстояниям (только для выбранных групп)
    multiple_min_stats = defaultdict(int)  # min_dist -> количество чисел с несколькими вариантами
    #Статистика по особым случаям
    special_cases_stats = defaultdict(int)  # min_dist -> количество особых случаев
    
    for chunk_result in results:
        for (i, code, min_dist, min_info_number, min_info_code, 
             all_min_codes, all_min_numbers, dist_distribution, all_distances) in chunk_result:
            if min_dist not in min_dist_counts:
                min_dist_counts[min_dist] = 0
            min_dist_counts[min_dist] += 1
            #Сохраняем распределение расстояний для этой группы
            for dist, count in dist_distribution.items():
                group_distributions[min_dist][dist] += count
                total_distances_per_group[min_dist] += count
            #Ищем особые случаи: расстояния больше текущего min_dist но меньше или равно заданному пользователем max_min_dist
            #ИСКЛЮЧАЕМ группу с min_dist равным заданному пользователем max_min_dist
            if min_dist in min_dist_groups and min_dist < user_max_min_dist:
                found_special_for_this_number = False
                #Проверяем все расстояния для данного числа
                for dist_idx, dist in enumerate(all_distances):
                    #Особый случай: расстояние больше текущего min_dist И меньше или равно заданному пользователем max_min_dist
                    if min_dist < dist <= user_max_min_dist:
                        special_info_number = inform_numbers[dist_idx]
                        special_info_code = inform_codes[dist_idx]
                        special_cases.append((
                            i, code, min_dist, min_info_number, min_info_code,
                            dist, special_info_number, special_info_code
                        ))
                        found_special_for_this_number = True
                #Если нашли особые случаи для этого числа, увеличиваем счетчик
                if found_special_for_this_number:
                    special_cases_stats[min_dist] += 1
            #Ищем случаи с несколькими кодами на минимальном расстоянии
            #ТОЛЬКО для выбранных групп min_dist
            if len(all_min_codes) > 1 and min_dist in min_dist_groups:
                multiple_min_cases.append((
                    i, code, min_dist, all_min_numbers, all_min_codes
                ))
                multiple_min_stats[min_dist] += 1
    
    execution_time = time.time() - start_time
    
    #Вывод основной статистики - показываем все группы от 1 до user_max_min_dist
    print(f"\nАнализ завершен за {execution_time:.2f} секунд")
    print(f"\nСтатистика минимальных расстояний для запрещенного диапазона:")
    print(f"(Показываются группы от 1 до {user_max_min_dist})")
    print("-" * 80)
    
    total_analyzed_numbers = sum(min_dist_counts.values())
    
    for dist in min_dist_groups:
        count = min_dist_counts.get(dist, 0)
        percentage = (count / total_analyzed_numbers) * 100 if total_analyzed_numbers > 0 else 0
        print(f"min_dist = {dist}: {count:,} чисел ({percentage:.2f}%)")
    
    #Детальный анализ распределения расстояний для выбранных групп
    print(f"\nДетальное распределение расстояний по группам:")
    print(f"(Показываются группы от 1 до {user_max_min_dist})")
    
    for min_dist_group in min_dist_groups:
        group_size = min_dist_counts.get(min_dist_group, 0)
        
        print(f"\nДля группы min_dist = {min_dist_group} ({group_size:,} чисел):")
        print("-" * 50)
        
        if group_size == 0:
            print("  Нет чисел в этой группе")
            continue
            
        print("Распределение расстояний до информационных кодов:")
        
        #Выводим все расстояния от 1 до total_mods
        for dist in range(1, total_mods + 1):
            count = group_distributions[min_dist_group].get(dist, 0)
            if total_distances_per_group[min_dist_group] > 0:
                percentage = (count / total_distances_per_group[min_dist_group]) * 100
                print(f"  Расстояние {dist}: {count:,} случаев ({percentage:.2f}%)")
            else:
                print(f"  Расстояние {dist}: {count:,} случаев (0.00%)")
        
        #Дополнительная статистика для группы
        total_in_group = total_distances_per_group[min_dist_group]
        if total_in_group > 0:
            avg_dist = sum(dist * count for dist, count in group_distributions[min_dist_group].items()) / total_in_group
            print(f"  Среднее расстояние в группе: {avg_dist:.2f}")
        print(f"  Всего сравнений в группе: {total_in_group:,}")
    
    #Вывод статистики по множественным минимальным расстояниям (только для выбранных групп)
    print(f"\nСтатистика по числам с несколькими информационными кодами на min_dist:")
    print(f"(Только для групп от 1 до {user_max_min_dist})")
    print("=" * 70)
    
    found_multiple_cases = False
    for min_dist_group in min_dist_groups:
        multiple_count = multiple_min_stats.get(min_dist_group, 0)
        group_size = min_dist_counts.get(min_dist_group, 0)
        
        if multiple_count > 0:
            found_multiple_cases = True
            percentage = (multiple_count / group_size) * 100 if group_size > 0 else 0
            print(f"min_dist = {min_dist_group}: {multiple_count:,} чисел ({percentage:.2f}% от группы)")
    
    if not found_multiple_cases:
        print("Не найдено чисел с несколькими информационными кодами на min_dist")
    
    #Вывод статистики по особым случаям
    print(f"\nСтатистика по особым случаям (расстояния > min_dist но <= {user_max_min_dist}):")
    print(f"(Только для групп от 1 до {user_max_min_dist-1})")
    print("=" * 70)
    
    total_special_cases_count = len(special_cases)
    found_special_cases = False
    
    for min_dist_group in min_dist_groups:
        if min_dist_group >= user_max_min_dist:  # Пропускаем группу с max_min_dist
            continue
            
        special_count = special_cases_stats.get(min_dist_group, 0)
        group_size = min_dist_counts.get(min_dist_group, 0)
        
        if special_count > 0:
            found_special_cases = True
            percentage = (special_count / group_size) * 100 if group_size > 0 else 0
            print(f"min_dist = {min_dist_group}: {special_count:,} чисел ({percentage:.2f}% от группы)")
    
    if not found_special_cases:
        print("Особые случаи не найдены")
    else:
        print(f"Всего особых случаев: {total_special_cases_count:,}")
    
    #Сохранение особых случаев
    if special_cases:
        print(f"\nНайдено особых случаев (расстояния > min_dist но <= {user_max_min_dist}): {len(special_cases):,}")
        save_choice = input('Сохранить эти особые случаи в файл? [yes/no]: ').strip().lower()
        
        if save_choice == 'yes':
            filename = 'special_cases_analysis.txt'
            save_special_cases(filename, special_cases, osn_inform)
            print(f"Особые случаи сохранены в файл: {filename}")
            
            #Показать несколько примеров на экране
            print(f"\nПримеры особых случаев (первые 5):")
            for i, case in enumerate(special_cases[:5]):
                (number, code, min_dist, min_info_number, min_info_code,
                 special_dist, special_info_number, special_info_code) = case
                
                print(f"Пример {i+1}:")
                print(f"Число из запрещенного диапазона: i = {number}, a = {code}")
                print(f"min_dist = {min_dist} (до i = {min_info_number}, a = {min_info_code})")
                print(f"Особое расстояние {special_dist} до i = {special_info_number}, a = {special_info_code}")
    else:
        print(f"\nОсобые случаи (расстояния > min_dist но <= {user_max_min_dist}) не найдены")
    
    #Сохранение случаев с несколькими кодами на min_dist
    if multiple_min_cases:
        print(f"\nНайдено случаев с несколькими кодами на min_dist: {len(multiple_min_cases):,}")
        save_choice = input('Сохранить случаи с несколькими кодами в файл? [yes/no]: ').strip().lower()
        
        if save_choice == 'yes':
            filename = 'multiple_min_cases_analysis.txt'
            save_multiple_min_cases(filename, multiple_min_cases, osn_inform)
            print(f"Случаи с несколькими кодами сохранены в файл: {filename}")
            
            #Показать несколько примеров на экране
            print(f"\nПримеры случаев с несколькими кодами на min_dist (первые 3):")
            print("=" * 70)
            for i, case in enumerate(multiple_min_cases[:3]):
                number, code, min_dist, all_min_numbers, all_min_codes = case
                
                print(f"Пример {i+1}:")
                print(f"Число из запрещенного диапазона: i = {number}, a = {code}")
                print(f"min_dist = {min_dist}")
                print(f"Количество информационных кодов на min_dist: {len(all_min_numbers)}")
                print(f"Информационные числа: {all_min_numbers}")
                print("-" * 40)
    else:
        print(f"\nСлучаи с несколькими кодами на min_dist не найдены")
    
    #Показать полную статистику по всем группам, если пользователь хочет
    show_full_stats = input('\nПоказать полную статистику по всем группам min_dist? [yes/no]: ').strip().lower()
    if show_full_stats == 'yes':
        print(f"\nПОЛНАЯ СТАТИСТИКА ПО ВСЕМ ГРУППАМ MIN_DIST:")
        print("=" * 50)
        total_all = sum(min_dist_counts.values())
        for dist in sorted(min_dist_counts.keys()):
            count = min_dist_counts[dist]
            percentage = (count / total_all) * 100 if total_all > 0 else 0
            
            #Показываем статистику по множественным кодам только для выбранных групп
            if dist in min_dist_groups:
                multiple_count = multiple_min_stats.get(dist, 0)
                multiple_percentage = (multiple_count / count) * 100 if count > 0 else 0
                
                #Показываем статистику по особым случаям только для групп меньше user_max_min_dist
                if dist < user_max_min_dist:
                    special_count = special_cases_stats.get(dist, 0)
                    special_percentage = (special_count / count) * 100 if count > 0 else 0
                    print(f"min_dist = {dist}: {count:,} чисел ({percentage:.2f}%)")
                    print(f"  - с несколькими кодами: {multiple_count:,} ({multiple_percentage:.2f}%)")
                    print(f"  - с особыми случаями: {special_count:,} ({special_percentage:.2f}%)")
                else:
                    print(f"min_dist = {dist}: {count:,} чисел ({percentage:.2f}%)")
                    print(f"  - с несколькими кодами: {multiple_count:,} ({multiple_percentage:.2f}%)")
            else:
                print(f"min_dist = {dist}: {count:,} чисел ({percentage:.2f}%)")

if __name__ == "__main__":
    analyze_modular_code_final()
