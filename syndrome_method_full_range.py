import numpy
import math
import time
import tracemalloc
from itertools import combinations, count, chain
from math import ceil
from collections import defaultdict

#Оценка размера файла
def estimate_file_size(full_diapazon, work_diapazon, osn):
    avg_line_length = len(f"i = 99999, a = {numpy.array([999]*len(osn))}, syndromes = {[999]*len(osn)}, corrected = {numpy.array([999]*len(osn))}, status = 'OK'\n")
    header_lines = 15
    avg_header_line = 50
    
    total_size = (full_diapazon * avg_line_length + 
                 header_lines * avg_header_line)
    
    return total_size

def convert_bytes(size):
    for x in ['bytes', 'KB', 'MB', 'GB']:
        if size < 1024.0:
            return "%3.1f %s" % (size, x)
        size /= 1024.0

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

def solve_linear_equation(equation_str, var_name, known_vars=None, mod=None):
    if known_vars is None:
        known_vars = {}
    
    for var, value in known_vars.items():
        equation_str = equation_str.replace(var, str(value))
    
    if mod is not None:
        for x in chain([0], (x*sign for x in count(1) for sign in (1, -1))):
            x_mod = x % mod
            expr = equation_str.replace(var_name, str(x_mod))
            try:
                value = eval(expr, {'__builtins__': None}, known_vars)
                if value % mod == 0:
                    return x_mod
            except:
                continue
    return None

#Вычисление синдромов для обнаруженных ошибок
def syndrome_decoding(osn_inform, osn_kontrol, chisl):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    k = len(osn_kontrol)
    n = len(osn_inform)
    n_plus_k = len(osn)
    B = calculate_orthogonal_bases(osn)
    B = numpy.array(B, dtype=object)
    
    B_modular = numpy.zeros((n_plus_k, n_plus_k), dtype=object)
    for i in range(n_plus_k):
        B_modular[i][:] = [b % int(osn[i]) for b in B]
    
    c = numpy.zeros((n_plus_k-1, n_plus_k), dtype=object)
    for i in range(n_plus_k-1):
        for j in range(n_plus_k-1):
            try:
                c[j][i] = pow(int(osn[i]), -1, int(osn[j+1]))
            except ValueError:
                continue
    
    for i in range(n_plus_k-2):
        c[i][i+1:] = 0
    
    B_opss = numpy.zeros((n_plus_k, n_plus_k), dtype=object)
    for k in range(n_plus_k):
        B_opss[k][0] = B_modular[k][0] % osn[0]
        if n_plus_k > 1:
            B_opss[k][1] = ((B_modular[k][1] - B_opss[k][0]) * c[0][0]) % osn[1]
        for i in range(2, n_plus_k):
            temp = (B_modular[k][i] - B_opss[k][0]) * c[i-1][0]
            for j in range(1, i):
                temp = (temp - B_opss[k][j]) * c[i-1][j]
            B_opss[k][i] = temp % osn[i]
    
    B_opss_chisl = numpy.zeros((n_plus_k, n_plus_k), dtype=object)
    for i in range(n_plus_k):
        for j in range(n_plus_k):
            B_opss_chisl[i][j] = (B_opss[i][j] * chisl[i]) if i < n else B_opss[i][j]
    
    opss = numpy.zeros(n_plus_k, dtype=object)
    for j in range(n_plus_k):
        for i in range(n_plus_k):
            if i < n:
                opss[j] += B_opss_chisl[i][j]
            else:
                if B_opss_chisl[i][j] != 0:
                    var_name = f"x{i-n+1}"
                    term = f"{B_opss_chisl[i][j]}*{var_name}"
                    opss[j] = f"{opss[j]} + {term}" if opss[j] != 0 else term
                    
    opss_total = numpy.zeros(n_plus_k, dtype=object)
    perenos = 0
    
    for i in range(n_plus_k):
        if isinstance(opss[i], (int, numpy.integer)):
            total = int(opss[i]) + int(perenos)
            opss_total[i] = total % osn[i]
            perenos = total // osn[i]
        else:
            opss_total[i] = f"({opss[i]} + {perenos})"
            perenos = f"({opss[i]} + {perenos}) // {osn[i]}"

    syndromes = []
    known_vars = {}
    
    for i in range(n, n_plus_k):
        mod = int(osn[i])
        equation = str(opss_total[i])
        var_name = f"x{i-n+1}"
        
        solution = solve_linear_equation(equation, var_name, known_vars, mod)
        
        if solution is not None:
            known_vars[var_name] = solution
            syndrome = (int(chisl[i]) - solution) % mod
            syndromes.append(syndrome)
        else:
            syndromes.append(None)

    return syndromes

#Исправление ошибок на основе синдромов
def locate_and_correct_errors(osn_inform, osn_kontrol, chisl, syndromes):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    k = len(osn_kontrol)
    n = len(osn_inform)
    work_diapazon = math.prod([int(x) for x in osn_inform])
    
    all_zero = all(s == 0 or s is None for s in syndromes)
    all_non_zero = all(s != 0 and s is not None for s in syndromes)
    some_non_zero = any(s != 0 and s is not None for s in syndromes) and not all_non_zero
    
    if all_zero:
        return chisl, "Нет ошибок", None
    
    elif all_non_zero:
        P = work_diapazon
        P_i = [P // int(x) for x in osn_inform]
        
        inv_P_i = []
        for i in range(n):
            inv_row = []
            for j in range(k):
                try:
                    inv = pow(int(P_i[i]), -1, int(osn[n+j]))
                    inv_row.append(inv)
                except ValueError:
                    inv_row.append(None)
            inv_P_i.append(inv_row)

        lambda_values_minus = []
        lambda_values_plus = []
        for i in range(n):
            lambda_row_minus = []
            lambda_row_plus = []
            for j in range(k):
                if inv_P_i[i][j] is not None and syndromes[j] is not None:
                    lambda_val_minus = (-syndromes[j] * inv_P_i[i][j]) % osn[n+j]
                    lambda_val_plus = (syndromes[j] * inv_P_i[i][j]) % osn[n+j]
                    lambda_row_minus.append(lambda_val_minus)
                    lambda_row_plus.append(lambda_val_plus)
                else:
                    lambda_row_minus.append(None)
                    lambda_row_plus.append(None)
            lambda_values_minus.append(lambda_row_minus)
            lambda_values_plus.append(lambda_row_plus)

        error_pos_minus = None
        error_value_minus = None
        error_pos_plus = None
        error_value_plus = None
        
        for i in range(n):
            unique_lambdas_minus = set(val for val in lambda_values_minus[i] if val is not None)
            if len(unique_lambdas_minus) == 1 and unique_lambdas_minus:
                error_pos_minus = i
                error_value_minus = (unique_lambdas_minus.pop() * P_i[i]) % osn_inform[i]
                break
            
        for i in range(n):
            unique_lambdas_plus = set(val for val in lambda_values_plus[i] if val is not None)
            if len(unique_lambdas_plus) == 1 and unique_lambdas_plus:
                error_pos_plus = i
                error_value_plus = (unique_lambdas_plus.pop() * P_i[i]) % osn_inform[i]
                break

        corrected_minus = None
        corrected_plus = None

        if error_pos_minus is not None:
            corrected_minus = numpy.array(chisl, dtype=object)
            corrected_minus[error_pos_minus] = (corrected_minus[error_pos_minus] - error_value_minus) % osn_inform[error_pos_minus]
        
        if error_pos_plus is not None:
            corrected_plus = numpy.array(chisl, dtype=object)
            corrected_plus[error_pos_plus] = (corrected_plus[error_pos_plus] + error_value_plus) % osn_inform[error_pos_plus]
            
        valid_corrections = []
        
        if corrected_minus is not None:
            A_minus = check_number(osn, corrected_minus)
            if A_minus < work_diapazon:
                valid_corrections.append(('minus', corrected_minus, error_value_minus, error_pos_minus))
        
        if corrected_plus is not None:
            A_plus = check_number(osn, corrected_plus)
            if A_plus < work_diapazon:
                valid_corrections.append(('plus', corrected_plus, error_value_plus, error_pos_plus))
        
        if valid_corrections:
            correction_type, corrected, error_value, error_pos = valid_corrections[0]
            message = f"Исправлено {'вычитанием' if correction_type == 'minus' else 'прибавлением'} {error_value} {'из' if correction_type == 'minus' else 'к'} позиции {error_pos}"
            
            #Вычисление вектора ошибки (всегда положительные значения)
            error_vector = numpy.zeros(len(chisl), dtype=object)
            if correction_type == 'minus':
                # Преобразуем отрицательную ошибку в положительную
                error_vector[error_pos] = (osn_inform[error_pos] - error_value) % osn_inform[error_pos]
            else:
                error_vector[error_pos] = error_value
            
            return corrected, message, error_vector
        else:
            return chisl, "Не удалось исправить ошибку", None
    
    elif some_non_zero:
        corrected = numpy.array(chisl, dtype=object)
        error_vector = numpy.zeros(len(chisl), dtype=object)
        for i in range(k):
            if syndromes[i] != 0 and syndromes[i] is not None:
                B = calculate_orthogonal_bases(osn)
                A = sum(int(chisl[j]) * int(B[j]) for j in range(n)) % work_diapazon
                old_value = corrected[n+i]
                corrected[n+i] = A % osn[n+i]
                error_vector[n+i] = (corrected[n+i] - old_value) % osn[n+i]
        return corrected, "Исправлена ошибка в контрольных основаниях", error_vector
    
    else:
        return chisl, "Не удалось определить тип ошибки", None

#Проверка методом проекций
def check_with_projections(osn_inform, osn_kontrol, chisl, osn):
    osn_kolvo = len(osn)
    work_diapazon = numpy.prod(osn_inform.astype(object))
    k = len(osn_kontrol)
    max_errors_to_correct = ceil(k / 2)
    
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
                    error_vector = numpy.zeros(len(chisl), dtype=object)
                    for i, pos in enumerate(error_positions):
                        old_value = corrected[pos]
                        corrected[pos] = solutions[i]
                        error_vector[pos] = (corrected[pos] - old_value) % osn_int[pos]
                    
                    A_corr = numpy.mod(numpy.sum(corrected * B_full), numpy.prod(osn_int.astype(object)))
                    if A_corr < work_diapazon:
                        return corrected, f"Исправлено {num_errors} ошибок методом проекций", error_vector
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
            return chisl, f"Обнаружено {num_errors} ошибок (исправление невозможно)", None
    
    return chisl, "Не удалось обнаружить ошибки (возможно, их больше чем k)", None

#Функция для вычисления разницы между векторами (всегда положительные значения)
def calculate_error_vector(original, corrected, osn):
    error_vector = numpy.zeros(len(original), dtype=object)
    for i in range(len(original)):
        # Всегда получаем положительное значение ошибки
        error_vector[i] = (corrected[i] - original[i]) % osn[i]
    return error_vector

#Функция для преобразования вектора ошибок в читаемый формат
def format_error_vector(error_vector):
    if error_vector is None:
        return "None"
    #Преобразуем np.int64 в обычные int для красивого вывода
    return tuple(int(x) if hasattr(x, 'item') else x for x in error_vector)

#Функция для преобразования массива в читаемый формат
def format_array(arr):
    if arr is None:
        return "None"
    # Преобразуем np.int64 в обычные int для красивого вывода
    return [int(x) if hasattr(x, 'item') else x for x in arr]

#Функция для сортировки векторов ошибок
def sort_error_vectors(error_vectors):
    """Сортирует векторы ошибок в логическом порядке: сначала по позиции, потом по величине ошибки"""
    return sorted(error_vectors, key=lambda vec: (
        #Сначала сортируем по позиции первой ненулевой ошибки
        next((i for i, x in enumerate(vec) if x != 0), len(vec)),
        #Затем по величине ошибки в этой позиции
        next((x for x in vec if x != 0), 0),
        #Затем по всему вектору
        vec
    ))

def main():
    # Начало программы
    start_time = time.time()
    tracemalloc.start()
    
    print('Основания нужно вводить в порядке возрастания: сначала информационные, затем контрольные. Они должны быть взаимно простыми')
    osn_inform_kolvo = int(input('Введите количество информационных оснований: '))
    osn_kontrol_kolvo = int(input('Введите количество контрольных оснований: '))
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

    #Оценка размера файла
    estimated_size = estimate_file_size(full_diapazon, work_diapazon, osn)
    print(f"\nПриблизительный размер файла с результатами: {convert_bytes(estimated_size)}")

    #Запрос на запись в файл
    filename = 'syndrome_decoding_results.txt'
    save_to_file = input(f'Записать все числа с проверками в файл {filename}? [yes/no]: ').strip().lower()

    file_created = False
    results = []
    syndrome_stats = defaultdict(list)
    error_vector_stats = defaultdict(list)

    if save_to_file == 'yes':
        file_created = True
        with open(filename, 'w', encoding='utf-8') as file:
            file.write(f'Информационные основания: {osn_inform}\n')
            file.write(f'Контрольные основания: {osn_kontrol}\n')
            file.write(f'Все основания: {osn}\n')
            file.write(f'Рабочий диапазон: 0-{work_diapazon-1}\n')
            file.write(f'Полный диапазон: 0-{full_diapazon-1}\n\n')
            
            numbers = generate_modular_numbers(osn_inform, osn_kontrol)
            
            file.write('Результаты синдромного декодирования\n')
            for i, a in numbers:
                syndromes = syndrome_decoding(osn_inform, osn_kontrol, a)
                corrected_syndrom, message_syndrom, error_vector_syndrom = locate_and_correct_errors(osn_inform, osn_kontrol, a, syndromes)
                corrected_proj, message_proj, error_vector_proj = check_with_projections(osn_inform, osn_kontrol, a, osn)
                
                #Если error_vector не был вычислен в функции, вычисляем его здесь
                if error_vector_syndrom is None:
                    error_vector_syndrom = calculate_error_vector(a, corrected_syndrom, osn)
                if error_vector_proj is None:
                    error_vector_proj = calculate_error_vector(a, corrected_proj, osn)
                
                #Сохраняем статистику по синдромам
                syndrome_tuple = tuple(syndromes)
                syndrome_stats[syndrome_tuple].append({
                    'i': i,
                    'original': a,
                    'corrected_syndrom': corrected_syndrom,
                    'error_vector_syndrom': error_vector_syndrom,
                    'corrected_proj': corrected_proj,
                    'error_vector_proj': error_vector_proj
                })
                
                #Сохраняем статистику по векторам ошибок
                error_vector_syndrom_tuple = tuple(format_error_vector(error_vector_syndrom))
                error_vector_stats[error_vector_syndrom_tuple].append({
                    'syndrome': syndromes,
                    'i': i,
                    'original': a,
                    'corrected': corrected_syndrom
                })
                
                #Проверка корректности
                A_syndrom = check_number(osn, corrected_syndrom)
                A_proj = check_number(osn, corrected_proj)
                
                status = 'OK' if A_syndrom == i and A_proj == i else 'ERROR'
                if A_syndrom > work_diapazon or A_proj > work_diapazon:
                    status = 'OUT_OF_RANGE'
                
                result_line = f"i = {i}, a = {format_array(a)}, syndromes = {syndromes}\n"
                result_line += f"  corrected_syndrom = {format_array(corrected_syndrom)}, error_vector_syndrom = {format_error_vector(error_vector_syndrom)}\n"
                result_line += f"  corrected_proj = {format_array(corrected_proj)}, error_vector_proj = {format_error_vector(error_vector_proj)}\n"
                result_line += f"  status = '{status}', message_syndrom = '{message_syndrom}', message_proj = '{message_proj}'\n\n"
                file.write(result_line)
                results.append((i, a, syndromes, corrected_syndrom, error_vector_syndrom, corrected_proj, error_vector_proj, status, message_syndrom, message_proj))
            
            #Записываем статистику по синдромам
            file.write('\n' + '='*80 + '\n')
            file.write('СТАТИСТИКА ПО СИНДРОМАМ И ВЕКТОРАМ ОШИБОК\n')
            file.write('='*80 + '\n\n')
            
            for syndrome, cases in syndrome_stats.items():
                file.write(f'Синдром: {syndrome}\n')
                file.write(f'Количество случаев: {len(cases)}\n')
                
                #Анализируем векторы ошибок
                error_vectors_syndrom = [case['error_vector_syndrom'] for case in cases]
                error_vectors_proj = [case['error_vector_proj'] for case in cases]
                
                #Находим уникальные векторы ошибок
                unique_errors_syndrom = set(tuple(format_error_vector(vec)) for vec in error_vectors_syndrom)
                unique_errors_proj = set(tuple(format_error_vector(vec)) for vec in error_vectors_proj)
                
                file.write(f'Уникальные векторы ошибок (синдромный метод): {len(unique_errors_syndrom)}\n')
                for error_vec in unique_errors_syndrom:
                    file.write(f'  {error_vec}\n')
                    # Находим пример для этого вектора ошибок
                    example_case = next(case for case in cases if tuple(format_error_vector(case['error_vector_syndrom'])) == error_vec)
                    file.write(f'    Пример: i={example_case["i"]}, a={format_array(example_case["original"])} -> corrected={format_array(example_case["corrected_syndrom"])}\n')
                
                file.write(f'Уникальные векторы ошибок (метод проекций): {len(unique_errors_proj)}\n')
                for error_vec in unique_errors_proj:
                    file.write(f'  {error_vec}\n')
                    # Находим пример для этого вектора ошибок
                    example_case = next(case for case in cases if tuple(format_error_vector(case['error_vector_proj'])) == error_vec)
                    file.write(f'    Пример: i={example_case["i"]}, a={format_array(example_case["original"])} -> corrected={format_array(example_case["corrected_proj"])}\n')
                
                file.write('\n')
            
            #Записываем статистику по векторам ошибок
            file.write('\n' + '='*80 + '\n')
            file.write('СТАТИСТИКА ПО ВЕКТОРАМ ОШИБОК И СИНДРОМАМ\n')
            file.write('='*80 + '\n\n')
            
            #Сортируем векторы ошибок
            sorted_error_vectors = sort_error_vectors([vec for vec in error_vector_stats.keys() if vec != tuple([0] * len(osn))])
            
            for error_vector in sorted_error_vectors:
                cases = error_vector_stats[error_vector]
                file.write(f'Вектор ошибок: {error_vector}\n')
                file.write(f'Количество случаев: {len(cases)}\n')
                
                #Находим уникальные синдромы для этого вектора ошибок
                unique_syndromes = set(tuple(case['syndrome']) for case in cases)
                file.write(f'Уникальные синдромы: {len(unique_syndromes)}\n')
                for syndrome in unique_syndromes:
                    file.write(f'  {syndrome}\n')
                    # Находим пример для этого синдрома
                    example_case = next(case for case in cases if tuple(case['syndrome']) == syndrome)
                    file.write(f'    Пример: i={example_case["i"]}, a={format_array(example_case["original"])} -> corrected={format_array(example_case["corrected"])}\n')
                
                file.write('\n')
        
        print(f'Файл успешно создан: {filename}')

    #Если не записываем в файл, все равно генерируем результаты
    if not file_created:
        numbers = generate_modular_numbers(osn_inform, osn_kontrol)
        for i, a in numbers:
            syndromes = syndrome_decoding(osn_inform, osn_kontrol, a)
            corrected_syndrom, message_syndrom, error_vector_syndrom = locate_and_correct_errors(osn_inform, osn_kontrol, a, syndromes)
            corrected_proj, message_proj, error_vector_proj = check_with_projections(osn_inform, osn_kontrol, a, osn)
            
            #Если error_vector не был вычислен в функции, вычисляем его здесь
            if error_vector_syndrom is None:
                error_vector_syndrom = calculate_error_vector(a, corrected_syndrom, osn)
            if error_vector_proj is None:
                error_vector_proj = calculate_error_vector(a, corrected_proj, osn)
            
            #Сохраняем статистику по синдромам
            syndrome_tuple = tuple(syndromes)
            syndrome_stats[syndrome_tuple].append({
                'i': i,
                'original': a,
                'corrected_syndrom': corrected_syndrom,
                'error_vector_syndrom': error_vector_syndrom,
                'corrected_proj': corrected_proj,
                'error_vector_proj': error_vector_proj
            })
            
            #Сохраняем статистику по векторам ошибок
            error_vector_syndrom_tuple = tuple(format_error_vector(error_vector_syndrom))
            error_vector_stats[error_vector_syndrom_tuple].append({
                'syndrome': syndromes,
                'i': i,
                'original': a,
                'corrected': corrected_syndrom
            })
            
            A_syndrom = check_number(osn, corrected_syndrom)
            A_proj = check_number(osn, corrected_proj)
            
            status = 'OK' if A_syndrom == i and A_proj == i else 'ERROR'
            if A_syndrom > work_diapazon or A_proj > work_diapazon:
                status = 'OUT_OF_RANGE'
            
            results.append((i, a, syndromes, corrected_syndrom, error_vector_syndrom, corrected_proj, error_vector_proj, status, message_syndrom, message_proj))

    #Запрос на вывод на экран
    user_input = input(f'\nВывести все числа с проверками на экран? (Их количество: {full_diapazon}) [yes/no]: ').strip().lower()

    if user_input == 'yes':
        print('\nРезультаты синдромного декодирования')
        for result in results:
            i, a, syndromes, corrected_syndrom, error_vector_syndrom, corrected_proj, error_vector_proj, status, message_syndrom, message_proj = result
            print(f"i = {i}, a = {format_array(a)}, syndromes = {syndromes}")
            print(f"  corrected_syndrom = {format_array(corrected_syndrom)}, error_vector_syndrom = {format_error_vector(error_vector_syndrom)}")
            print(f"  corrected_proj = {format_array(corrected_proj)}, error_vector_proj = {format_error_vector(error_vector_proj)}")
            print(f"  status = '{status}', message_syndrom = '{message_syndrom}', message_proj = '{message_proj}'")
            print()

    #Вывод статистики по синдромам
    print('\n' + '='*80)
    print('СТАТИСТИКА ПО СИНДРОМАМ И ВЕКТОРАМ ОШИБОК')
    print('='*80)
    
    for syndrome, cases in syndrome_stats.items():
        print(f'\nСиндром: {syndrome}')
        print(f'Количество случаев: {len(cases)}')
        
        #Анализируем векторы ошибок
        error_vectors_syndrom = [case['error_vector_syndrom'] for case in cases]
        error_vectors_proj = [case['error_vector_proj'] for case in cases]
        
        #Находим уникальные векторы ошибок
        unique_errors_syndrom = set(tuple(format_error_vector(vec)) for vec in error_vectors_syndrom)
        unique_errors_proj = set(tuple(format_error_vector(vec)) for vec in error_vectors_proj)
        
        print(f'Уникальные векторы ошибок (синдромный метод): {len(unique_errors_syndrom)}')
        for error_vec in unique_errors_syndrom:
            print(f'  {error_vec}')
            # Находим пример для этого вектора ошибок
            example_case = next(case for case in cases if tuple(format_error_vector(case['error_vector_syndrom'])) == error_vec)
            print(f'    Пример: i={example_case["i"]}, a={format_array(example_case["original"])} -> corrected={format_array(example_case["corrected_syndrom"])}')
        
        print(f'Уникальные векторы ошибок (метод проекций): {len(unique_errors_proj)}')
        for error_vec in unique_errors_proj:
            print(f'  {error_vec}')
            #Находим пример для этого вектора ошибок
            example_case = next(case for case in cases if tuple(format_error_vector(case['error_vector_proj'])) == error_vec)
            print(f'    Пример: i={example_case["i"]}, a={format_array(example_case["original"])} -> corrected={format_array(example_case["corrected_proj"])}')

    #Вывод статистики по векторам ошибок
    print('\n' + '='*80)
    print('СТАТИСТИКА ПО ВЕКТОРАМ ОШИБОК И СИНДРОМАМ')
    print('='*80)
    
    #Сортируем векторы ошибок
    sorted_error_vectors = sort_error_vectors([vec for vec in error_vector_stats.keys() if vec != tuple([0] * len(osn))])
    
    for error_vector in sorted_error_vectors:
        cases = error_vector_stats[error_vector]
        print(f'\nВектор ошибок: {error_vector}')
        print(f'Количество случаев: {len(cases)}')
        
        #Находим уникальные синдромы для этого вектора ошибок
        unique_syndromes = set(tuple(case['syndrome']) for case in cases)
        print(f'Уникальные синдромы: {len(unique_syndromes)}')
        for syndrome in unique_syndromes:
            print(f'  {syndrome}')
            # Находим пример для этого синдрома
            example_case = next(case for case in cases if tuple(case['syndrome']) == syndrome)
            print(f'    Пример: i={example_case["i"]}, a={format_array(example_case["original"])} -> corrected={format_array(example_case["corrected"])}')

    #Замер памяти и времени
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    execution_time = time.time() - start_time

    print('\nРесурсы программы')
    print(f'Использовано памяти: {current/1024:.2f} KB')
    print(f'Пиковое использование памяти: {peak/1024:.2f} KB')
    print(f'Время выполнения: {execution_time:.2f} секунд')

    if file_created:
        print(f'\nПрограмма завершена. Результаты сохранены в файле: {filename}')
    else:
        print('\nПрограмма завершена. Файл не создавался.')

if __name__ == "__main__":
    main()
