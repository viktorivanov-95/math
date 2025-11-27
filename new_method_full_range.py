import numpy
import math
import time
import tracemalloc
from itertools import combinations, count, chain
from math import ceil
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import statistics
import pandas as pd
import sys
from tqdm import tqdm

#Генерация модулярных чисел
def generate_modular_numbers(osn_inform, osn_kontrol):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    full_diapazon = numpy.prod(osn.astype(object))
    numbers = []
    for i in range(int(full_diapazon)):
        a = numpy.mod(i, osn)
        numbers.append((i, a))
    return numbers

#Вычисление ортогональных базисов
def calculate_orthogonal_bases(osn):
    full_diapazon = math.prod([int(x) for x in osn])
    P = [full_diapazon // int(x) for x in osn]
    beta = [p % m for p, m in zip(P, osn)]
    m = [pow(int(b), -1, int(m)) if b != 0 else 0 for b, m in zip(beta, osn)]
    return [mi * pi for mi, pi in zip(m, P)]

#Проверка числа методом ортогональных базисов
def check_number(osn, number):
    B = calculate_orthogonal_bases(osn)
    total = sum(int(n) * int(b) for n, b in zip(number, B))
    A = total % math.prod([int(x) for x in osn])
    return A

#Вычисление обнаружения ошибки новым методом:
def new_method_decoding(osn_inform, osn_kontrol, chisl, max_error_correction):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    k = len(osn_kontrol)
    n = len(osn_inform)
    n_plus_k = len(osn)
    full_diapazon = math.prod([int(x) for x in osn])
    osn_inform_mult = math.prod([int(x) for x in osn_inform])
    osn_kontrol_mult = math.prod([int(x) for x in osn_kontrol])
    A = check_number(osn, chisl)
    if A < osn_inform_mult:
        return chisl.tolist()
    else:
        #Проверка на ошибки только по информационным модулям
        inform_errors = A - int((A/osn_kontrol_mult))*osn_kontrol_mult
        inform_error_vector = [(int((A/osn_kontrol_mult))*osn_kontrol_mult) % osn[i] for i in range(len(osn))]
        inform_error_vector_checking = sum(1 for x in inform_error_vector if x != 0)
        if inform_errors < osn_inform_mult and inform_error_vector_checking <= max_error_correction:
            #Вычитаем из самого числа вектор ошибки:
            right_modular_vector = [(chisl[i] - inform_error_vector[i]) % osn[i] for i in range(len(osn))]
            A = check_number(osn, right_modular_vector)
            if A < osn_inform_mult:
                return right_modular_vector
            else:
                return None
        elif inform_errors >= osn_inform_mult or inform_error_vector_checking > max_error_correction:
            #Проверка на ошибки только по контрольным модулям
            control_error_vector = [(int((A/osn_inform_mult))*osn_inform_mult) % osn[i] for i in range(len(osn))]
            control_error_vector_checking = sum(1 for x in control_error_vector if x != 0)
            if control_error_vector_checking <= max_error_correction:
                right_modular_vector = [(chisl[i] - control_error_vector[i]) % osn[i] for i in range(len(osn))]
                A = check_number(osn, right_modular_vector)
                if A < osn_inform_mult:
                    return right_modular_vector
                else:
                    return None
            elif control_error_vector_checking > max_error_correction and max_error_correction !=1:
                #Проверка на ошибки по информационным и контрольным модулям. max_error_correction должен быть больше 1
                #Вычисляем наборы контрольных оснований для проверок
                # Вычисляем диапазон размеров комбинаций
                max_size = len(osn_kontrol) - (len(osn_kontrol) - max_error_correction)
                min_size = max_error_correction-1
                if min_size == 1:
                    min_size = min_size + 1
                control_modules_for_checking = []
                for size in range(min_size, max_size + 1):
                    for combination in combinations(osn_kontrol, size):
                        control_osn_combiation_mult = math.prod(combination)
                        control_modules_for_checking.append(control_osn_combiation_mult)
                for combination in control_modules_for_checking:
                    combined_error_vector = [(int((A/combination))*combination) % osn[i] for i in range(len(osn))]
                    combined_error_vector_checking = sum(1 for x in combined_error_vector if x != 0)
                    #Следующие два условия для случая, если произведение первых двух контрольных оснований приближенно равны
                    #произведению информационных оснований
                    if combined_error_vector == chisl.tolist() and combined_error_vector_checking > max_error_correction:
                        combined_error_vector = [((int(A/combination)*(combination)) - combination) % osn[i] for i in range(len(osn))]
                        combined_error_vector_checking = sum(1 for x in combined_error_vector if x != 0)
                    elif combined_error_vector == chisl.tolist() and combined_error_vector_checking <= max_error_correction:
                        combined_error_vector = [((int(A/combination)*(combination))) % osn[i] for i in range(len(osn))]
                        combined_error_vector_checking = sum(1 for x in combined_error_vector if x != 0)
                    if combined_error_vector_checking <= max_error_correction and (A - int(A/combination)*combination) < osn_inform_mult:
                        right_modular_vector = [(chisl[i] - combined_error_vector[i]) % osn[i] for i in range(len(osn))]
                        check_number_orth_bas = check_number(osn, right_modular_vector)
                        if check_number_orth_bas < osn_inform_mult:
                            return right_modular_vector
                        else:
                            continue
                    if combined_error_vector_checking > max_error_correction and control_modules_for_checking.index(combination) == k-1:
                        return None
            elif control_error_vector_checking > max_error_correction and max_error_correction !=1:
                return None
    return None


#Метод проекций
def check_with_projections(osn_inform, osn_kontrol, chisl, osn, max_error_correction):
    osn_kolvo = len(osn)
    work_diapazon = numpy.prod(osn_inform.astype(object))
    k = len(osn_kontrol)
    max_errors_to_correct = max_error_correction
    osn_int = osn.astype(int)
    for num_errors in range(1, max_errors_to_correct + 1):
        for error_positions in combinations(range(osn_kolvo), num_errors):
            proj_osn = numpy.delete(osn_int, error_positions)
            proj_chisl = numpy.delete(chisl, error_positions)
            B_proj = calculate_orthogonal_bases(proj_osn)
            A_proj = numpy.mod(numpy.sum(proj_chisl * B_proj), numpy.prod(proj_osn.astype(object)))
            if A_proj < work_diapazon:
                B_full = calculate_orthogonal_bases(osn_int)
                M = numpy.zeros((num_errors, num_errors), dtype=float)
                for i, pos in enumerate(error_positions):
                    for j, p in enumerate(error_positions):
                        M[i,j] = float(B_full[pos] % osn_int[p])
                b = numpy.zeros(num_errors, dtype=float)
                for i, pos in enumerate(error_positions):
                    s = sum(chisl[p] * B_full[p] for p in range(osn_kolvo) if p not in error_positions)
                    b[i] = float((A_proj - s) % osn_int[pos])
                try:
                    solutions = numpy.linalg.solve(M, b)
                    solutions = numpy.mod(numpy.round(solutions).astype(int), osn_int[list(error_positions)])
                    corrected = chisl.copy()
                    for i, pos in enumerate(error_positions):
                        corrected[pos] = int(solutions[i])
                    A_corr = numpy.mod(numpy.sum(corrected * B_full), numpy.prod(osn_int.astype(object)))
                    if A_corr < work_diapazon:
                        return corrected.tolist()
                except numpy.linalg.LinAlgError:
                    continue
    
    for num_errors in range(1, k + 1):
        all_valid = True
        for error_positions in combinations(range(osn_kolvo), num_errors):
            temp_chisl = chisl.copy()
            for pos in error_positions:
                temp_chisl[pos] = (temp_chisl[pos] + 1) % osn_int[pos]
                
            B_full = calculate_orthogonal_bases(osn_int)
            A_err = numpy.mod(numpy.sum(temp_chisl * B_full), numpy.prod(osn_int.astype(object)))
            if A_err < work_diapazon:
                all_valid = False
                break
        
        if all_valid:
            return None
    
    return None

#Функция для обработки одного числа
def process_single_number(args):
    i, a, osn_inform, osn_kontrol, osn, work_diapazon, max_error_correction = args
    #Применяем оба метода
    new_method_result = new_method_decoding(osn_inform, osn_kontrol, a, max_error_correction)
    projection_result = check_with_projections(osn_inform, osn_kontrol, a, osn, max_error_correction)
    #Проверяем корректность результатов
    A_new = check_number(osn, new_method_result) if new_method_result is not None else None
    A_proj = check_number(osn, projection_result) if projection_result is not None else None
    #Определяем статус
    new_method_ok = (A_new == i) if new_method_result is not None else False
    projection_ok = (A_proj == i) if projection_result is not None else False
    #Проверяем, находятся ли результаты в рабочем диапазоне
    new_method_in_range = (A_new < work_diapazon) if new_method_result is not None else False
    projection_in_range = (A_proj < work_diapazon) if projection_result is not None else False
    #Проверяем совпадение результатов (оба None считаются совпадением)
    if new_method_result is None and projection_result is None:
        results_match = True
    elif new_method_result is not None and projection_result is not None:
        results_match = (new_method_result == projection_result)
    else:
        results_match = False
    return {
        'i': i,
        'a': a.tolist(),
        'new_method_result': new_method_result,
        'projection_result': projection_result,
        'new_method_ok': new_method_ok,
        'projection_ok': projection_ok,
        'new_method_in_range': new_method_in_range,
        'projection_in_range': projection_in_range,
        'results_match': results_match,
        'A_new': A_new,
        'A_proj': A_proj
    }

#Оценка размера файла
def estimate_file_size(num_discrepancies, avg_line_length=200):
    return num_discrepancies * avg_line_length

def convert_bytes(size):
    for x in ['bytes', 'KB', 'MB', 'GB']:
        if size < 1024.0:
            return "%3.1f %s" % (size, x)
        size /= 1024.0

#Сравнение алгоритмов
def scientific_comparison(results, execution_time, max_error_correction):
    
    print('\n' + '='*100)
    print('Сравнение алгоритмов декодирования')
    print('='*100)
    print(f'Максимальная кратность ошибки: {max_error_correction}')
    print('='*100)
    
    #Базовые метрики
    total_numbers = len(results)
    new_method_correct = sum(1 for r in results if r['new_method_ok'])
    projection_correct = sum(1 for r in results if r['projection_ok'])
    new_method_in_range = sum(1 for r in results if r['new_method_in_range'])
    projection_in_range = sum(1 for r in results if r['projection_in_range'])
    results_match = sum(1 for r in results if r['results_match'])
    discrepancies = [r for r in results if not r['results_match']]
    
    #Эффективность исправления ошибок
    print('\n1. Эффективность исправления ошибок:')
    print('-' * 50)
    
    new_method_accuracy = new_method_correct / total_numbers * 100
    projection_accuracy = projection_correct / total_numbers * 100
    
    print(f'Точность нового метода: {new_method_accuracy:.4f}% ({new_method_correct}/{total_numbers})')
    print(f'Точность метода проекций: {projection_accuracy:.4f}% ({projection_correct}/{total_numbers})')
    
    if new_method_accuracy > projection_accuracy:
        advantage = new_method_accuracy - projection_accuracy
        print(f'Новый метод превосходит на: {advantage:.4f}%')
    elif projection_accuracy > new_method_accuracy:
        advantage = projection_accuracy - new_method_accuracy
        print(f'Метод проекций превосходит на: {advantage:.4f}%')
    else:
        print('Методы демонстрируют одинаковую точность')
    
    #Анализ расхождений
    print('\n2. Анализ расхождений:')
    print('-' * 50)
    
    if discrepancies:
        #Анализ типов расхождений
        new_method_better = 0
        projection_better = 0
        both_wrong = 0
        different_corrections = 0
        one_none_other_not = 0
        
        for disc in discrepancies:
            new_ok = disc['new_method_ok']
            proj_ok = disc['projection_ok']
            new_result = disc['new_method_result']
            proj_result = disc['projection_result']
            
            #Оба метода вернули None - это не расхождение (уже отфильтровано)
            if new_result is None and proj_result is None:
                continue
                
            #Один метод вернул None, другой - значение
            elif new_result is None or proj_result is None:
                one_none_other_not += 1
                if new_ok and not proj_ok:
                    new_method_better += 1
                elif proj_ok and not new_ok:
                    projection_better += 1
                elif not new_ok and not proj_ok:
                    both_wrong += 1
                    
            #Оба метода вернули значения, но разные
            elif new_ok and not proj_ok:
                new_method_better += 1
            elif proj_ok and not new_ok:
                projection_better += 1
            elif not new_ok and not proj_ok:
                both_wrong += 1
            else:  # оба правильные, но разные результаты
                different_corrections += 1
        
        print(f'Случаев, когда новый метод лучше: {new_method_better}')
        print(f'Случаев, когда метод проекций лучше: {projection_better}')
        print(f'Случаев, когда оба метода ошибаются: {both_wrong}')
        print(f'Случаев разных корректных исправлений: {different_corrections}')
        print(f'Случаев, когда один метод вернул None, а другой - значение: {one_none_other_not}')
        
        #Эффективность при наличии расхождений
        if new_method_better + projection_better > 0:
            new_method_superiority = new_method_better / (new_method_better + projection_better) * 100
            print(f'Эффективность нового метода в спорных случаях: {new_method_superiority:.2f}%')
    else:
        print('Расхождения между методами отсутствуют')
    
    #Производительность
    print('\n3. Производительность:')
    print('-' * 50)
    print(f'Общее время выполнения: {execution_time:.4f} секунд')
    
    return {
        'new_method_accuracy': new_method_accuracy,
        'projection_accuracy': projection_accuracy,
        'total_discrepancies': len(discrepancies),
        'execution_time': execution_time,
        'max_error_correction': max_error_correction
    }

#Оптимальное количество рабочих процессов
def calculate_optimal_workers():
    total_cores = multiprocessing.cpu_count()
    optimal_workers = max(1, int(total_cores * 0.75))
    return optimal_workers

#Обрабатывает пакет чисел
def process_batch(batch_args):
    return [process_single_number(args) for args in batch_args]

def main():
    start_time = time.time()
    tracemalloc.start()
    
    print('\nОснования нужно вводить в порядке возрастания: сначала информационные, затем контрольные. Они должны быть взаимно простыми')
    osn_inform_kolvo = int(input('Введите количество информационных оснований: '))
    osn_kontrol_kolvo = int(input('Введите количество контрольных оснований: '))
    max_error_correction = int(input('Введите максимальную кратность ошибки: '))
    osn_inform = numpy.zeros((osn_inform_kolvo), dtype=object)
    osn_kontrol = numpy.zeros((osn_kontrol_kolvo), dtype=object)
    for i in range(osn_inform_kolvo):
        osn_inform[i] = int(input(f'Введите {i+1} информационное основание: '))
    for i in range(osn_kontrol_kolvo):
        osn_kontrol[i] = int(input(f'Введите {i+1} контрольное основание: '))

    osn = numpy.concatenate((osn_inform, osn_kontrol))
    work_diapazon = numpy.prod(osn_inform.astype(object))
    full_diapazon = numpy.prod(osn.astype(object))

    print(f'\nРабочий диапазон: 0-{work_diapazon-1} ({work_diapazon} чисел)')
    print(f'Полный диапазон: 0-{full_diapazon-1} ({full_diapazon} чисел)')
    print(f'Максимальная кратность ошибки: {max_error_correction}')
    
    #Информация о системе
    total_cores = multiprocessing.cpu_count()
    optimal_workers = calculate_optimal_workers()
    print(f'\nИнформация о системе:')
    print(f'Всего ядер процессора: {total_cores}')
    print(f'Рекомендуемое количество рабочих процессов: {optimal_workers}')
    
    #Запрос на использование параллельной обработки
    use_parallel = input('Использовать параллельную обработку? [yes/no]: ').strip().lower() == 'yes'
    
    if use_parallel:
        custom_workers = input(f'Использовать рекомендуемое количество процессов ({optimal_workers}) или указать свое? [recommend/custom]: ').strip().lower()
        if custom_workers == 'custom':
            optimal_workers = int(input(f'Введите количество процессов (1-{total_cores}): '))
            optimal_workers = max(1, min(optimal_workers, total_cores))
    
    #Генерация всех чисел
    print(f'\nГенерация модулярных чисел...')
    numbers = generate_modular_numbers(osn_inform, osn_kontrol)
    print(f'Сгенерировано {len(numbers)} чисел для проверки')
    
    #Подготовка аргументов для обработки
    args_list = [(i, a, osn_inform, osn_kontrol, osn, work_diapazon, max_error_correction) for i, a in numbers]
    
    results = []
    
    print(f'\nНачало обработки...')
    processing_start = time.time()
    
    if use_parallel:
        print(f'Используется параллельная обработка ({optimal_workers} процессов)')
        
        #Разбиваем на пакеты для лучшего отображения прогресса
        batch_size = max(1, len(args_list) // 100)  # 100 пакетов или меньше
        batches = [args_list[i:i + batch_size] for i in range(0, len(args_list), batch_size)]
        
        with ProcessPoolExecutor(max_workers=optimal_workers) as executor:
            #Отправляем все пакеты на выполнение
            future_to_batch = {executor.submit(process_batch, batch): i for i, batch in enumerate(batches)}
            
            #Создаем прогресс-бар
            with tqdm(total=len(batches), desc="Обработка пакетов", unit="пакет") as pbar:
                for future in as_completed(future_to_batch):
                    batch_results = future.result()
                    results.extend(batch_results)
                    pbar.update(1)
                    pbar.set_postfix({
                        'обработано': f"{len(results)}/{len(args_list)}",
                        'время': f"{time.time() - processing_start:.1f}с"
                    })
                    
    else:
        print('Используется последовательная обработка')
        #Создаем прогресс-бар для последовательной обработки
        for args in tqdm(args_list, desc="Обработка чисел", unit="число"):
            result = process_single_number(args)
            results.append(result)
    
    processing_time = time.time() - processing_start
    print(f'\nОбработка завершена за {processing_time:.2f} секунд')
    
    #Общее время выполнения
    execution_time = time.time() - start_time
    
    #Научное сравнение
    comparison_results = scientific_comparison(results, execution_time, max_error_correction)
    
    #Базовая статистика для быстрого обзора
    discrepancies = [r for r in results if not r['results_match']]
    real_discrepancies = [r for r in discrepancies if not (r['new_method_result'] is None and r['projection_result'] is None)]
    
    print('\n' + '='*80)
    print('Базовая статистика результатов')
    print('='*80)
    print(f'Всего проверено чисел: {len(results)}')
    print(f'Корректно исправлено новым методом: {sum(1 for r in results if r["new_method_ok"])} ({comparison_results["new_method_accuracy"]:.2f}%)')
    print(f'Корректно исправлено методом проекций: {sum(1 for r in results if r["projection_ok"])} ({comparison_results["projection_accuracy"]:.2f}%)')
    print(f'Реальных расхождений (исключая оба None): {len(real_discrepancies)} ({len(real_discrepancies)/len(results)*100:.2f}%)')
    
    #Вывод первых 5 реальных расхождений
    if real_discrepancies:
        print(f'\nПервые 5 расхождений:')
        for i, disc in enumerate(real_discrepancies[:5]):
            print(f'{i+1}. i = {disc["i"]}, a = {disc["a"]}')
            print(f'   Новый метод: {disc["new_method_result"]} (A = {disc["A_new"]}, OK: {disc["new_method_ok"]})')
            print(f'   Метод проекций: {disc["projection_result"]} (A = {disc["A_proj"]}, OK: {disc["projection_ok"]})')
            print()
    else:
        print('\nРеальных расхождений не обнаружено')
    
    #Предложение сохранить все реальные расхождения в файл
    if real_discrepancies:
        estimated_size = estimate_file_size(len(real_discrepancies))
        print(f'Приблизительный размер файла со всеми РЕАЛЬНЫМИ расхождениями: {convert_bytes(estimated_size)}')
        
        save_discrepancies = input('Сохранить все РЕАЛЬНЫЕ расхождения в файл? [yes/no]: ').strip().lower()
        
        if save_discrepancies == 'yes':
            filename = 'real_discrepancies_analysis.txt'
            print(f'Сохранение в файл {filename}')
            
            with open(filename, 'w', encoding='utf-8') as file:
                file.write('Анализ расхождений между методами декодирования\n')
                
                file.write(f'Информационные основания: {osn_inform.tolist()}\n')
                file.write(f'Контрольные основания: {osn_kontrol.tolist()}\n')
                file.write(f'Все основания: {osn.tolist()}\n')
                file.write(f'Рабочий диапазон: 0-{work_diapazon-1}\n')
                file.write(f'Полный диапазон: 0-{full_diapazon-1}\n')
                file.write(f'Максимальная кратность ошибки: {max_error_correction}\n\n')
                
                file.write('Статистика сравнения:\n')
                file.write(f'Точность нового метода: {comparison_results["new_method_accuracy"]:.4f}%\n')
                file.write(f'Точность метода проекций: {comparison_results["projection_accuracy"]:.4f}%\n')
                file.write(f'Реальных расхождений: {len(real_discrepancies)}\n')
                file.write(f'Время выполнения: {execution_time:.4f} секунд\n\n')
                
                file.write('ДЕТАЛЬНЫЙ АНАЛИЗ РЕАЛЬНЫХ РАСХОЖДЕНИЙ:\n')
                file.write('='*80 + '\n\n')
                
                for i, disc in enumerate(real_discrepancies):
                    file.write(f'Расхождение {i+1}:\n')
                    file.write(f'Исходное число: i = {disc["i"]}, a = {disc["a"]}\n')
                    file.write(f'Новый метод: {disc["new_method_result"]} (A = {disc["A_new"]}, OK: {disc["new_method_ok"]})\n')
                    file.write(f'Метод проекций:n {disc["projection_result"]} (A = {disc["A_proj"]}, OK: {disc["projection_ok"]})\n')
                    
                    # Анализ конкретного случая
                    if disc['new_method_ok'] and not disc['projection_ok']:
                        file.write('Новый метод корректен, метод проекций ошибается\n')
                    elif disc['projection_ok'] and not disc['new_method_ok']:
                        file.write('Метод проекций корректен, новый метод ошибается\n')
                    elif not disc['new_method_ok'] and not disc['projection_ok']:
                        file.write('Оба метода ошибаются\n')
                    else:
                        file.write('Оба метода корректны, но дают разные результаты\n')
                    
                    file.write('\n')
            
            print(f'Подробный анализ расхождений сохранен в файл: {filename}')
    
    #Замер памяти
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    print('\n Ресурсы программы:')
    print(f'Использовано памяти: {current/1024:.2f} KB')
    print(f'Пиковое использование памяти: {peak/1024:.2f} KB')
    print(f'Общее время выполнения: {execution_time:.2f} секунд')
    print(f'Чистое время обработки: {processing_time:.2f} секунд')

if __name__ == "__main__":
    main()
