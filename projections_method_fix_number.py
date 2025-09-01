import math
from itertools import combinations
from math import ceil, prod
import time
import tracemalloc

#Вычисление ортогональных базисов
def calculate_orthogonal_bases(osn):
    print('Вычисление ортогональных базисов')
    full_diapazon = prod(osn)
    print(f"1. Полный диапазон P = {full_diapazon}")

    P = [full_diapazon // p for p in osn]
    print(f"2. Величины P/pi:")
    for i in range(len(osn)):
        print(f"P[{i}] = {full_diapazon} // {osn[i]} = {P[i]}")
    
    beta = [P_i % p for P_i, p in zip(P, osn)]
    print(f"3. Величины βi = Pi mod pi:")
    for i in range(len(osn)):
        print(f"β[{i}] = {P[i]} mod {osn[i]} = {beta[i]}")
    
    print(f"4. Обратные элементы mi = inv(βi) mod pi:")
    m = []
    for i in range(len(osn)):
        try:
            m_val = pow(int(beta[i]), -1, int(osn[i]))
            m.append(m_val)
            print(f"m[{i}] = inv({beta[i]}) mod {osn[i]} = {m_val}")
        except ValueError:
            print(f"Ошибка: невозможно найти обратный элемент для β[{i}] = {beta[i]} mod {osn[i]}")
            return None
    
    result = [m_i * P_i for m_i, P_i in zip(m, P)]
    print(f"5. Ортогональные базисы Bi = mi * Pi:")
    for i in range(len(osn)):
        print(f"B[{i}] = {m[i]} * {P[i]} = {result[i]}")
    
    return result

# Безопасное вычисление суммы с модульной арифметикой
def safe_modular_sum(values, mod):
    """Вычисляет сумму значений по модулю, избегая переполнения"""
    result = 0
    for value in values:
        result = (result + value) % mod
    return result

# Безопасное вычисление скалярного произведения с модульной арифметикой
def safe_dot_product(a, b, mod):
    """Вычисляет скалярное произведение по модулю, избегая переполнения"""
    result = 0
    for a_i, b_i in zip(a, b):
        # Вычисляем каждое слагаемое по модулю
        term = (a_i % mod) * (b_i % mod) % mod
        result = (result + term) % mod
    return result

# Решение системы уравнений в модульной арифметике
def solve_modular_system(M, b, mods):
    n = len(b)
    print(f"Решение системы {n}x{n} в модульной арифметике")
    
    # Для одиночной ошибки
    if n == 1:
        print(f"Одиночная ошибка: M[0,0] = {M[0][0]}, b[0] = {b[0]}, mod = {mods[0]}")
        try:
            inv = pow(int(M[0][0]), -1, int(mods[0]))
            print(f"Обратный элемент: inv({M[0][0]}) mod {mods[0]} = {inv}")
            solution = (b[0] * inv) % mods[0]
            print(f"Решение: ({b[0]} * {inv}) mod {mods[0]} = {solution}")
            return [solution]
        except Exception as e:
            print(f"Ошибка при вычислении обратного элемента: {e}")
            return None
    
    # Для множественных ошибок используем метод Гаусса с модульной арифметикой
    print("Метод Гаусса для множественных ошибок")
    M = [row[:] for row in M]  # Копируем матрицу
    b = b[:]  # Копируем вектор
    
    for i in range(n):
        print(f"\nШаг {i+1}:")
        print(f"Матрица M:")
        for row in M:
            print(row)
        print(f"Вектор b: {b}")
        
        # Находим строку с ненулевым ведущим элементом
        if M[i][i] == 0:
            print(f"Ведущий элемент M[{i}][{i}] = 0, ищем замену...")
            for j in range(i + 1, n):
                if M[j][i] != 0:
                    print(f"Меняем строки {i} и {j}")
                    # Меняем строки местами
                    M[i], M[j] = M[j], M[i]
                    b[i], b[j] = b[j], b[i]
                    print(f"После замены:")
                    for row in M:
                        print(row)
                    print(f"b: {b}")
                    break
        
        if M[i][i] == 0:
            print("Система вырождена")
            return None
        
        print(f"Ведущий элемент M[{i}][{i}] = {M[i][i]}")
        
        # Нормализуем текущую строку
        try:
            inv = pow(int(M[i][i]), -1, int(mods[i]))
            print(f"Обратный элемент: inv({M[i][i]}) mod {mods[i]} = {inv}")
            
            # Нормализуем строку матрицы
            for j in range(n):
                M[i][j] = (M[i][j] * inv) % mods[i]
            b[i] = (b[i] * inv) % mods[i]
            
            print(f"После нормализации:")
            for row in M:
                print(row)
            print(f"b: {b}")
        except Exception as e:
            print(f"Ошибка при нормализации: {e}")
            return None
        
        # Исключаем текущую переменную из других уравнений
        for j in range(n):
            if j != i:
                factor = M[j][i]
                print(f"Исключаем из строки {j}, factor = {factor}")
                
                # Вычитаем из строки j
                for k in range(n):
                    M[j][k] = (M[j][k] - factor * M[i][k]) % mods[j]
                b[j] = (b[j] - factor * b[i]) % mods[j]
                
                print(f"После исключения:")
                for row in M:
                    print(row)
                print(f"b: {b}")
    
    print(f"Конечный результат:")
    for row in M:
        print(row)
    print(f"b: {b}")
    return b

#Вычисление методом проекций
def check_with_projections(osn_inform, osn_kontrol, chisl, osn, max_errors_mode):
    print('Вычисление методом проекций')
    osn_kolvo = len(osn)
    work_diapazon = prod(osn_inform)
    k = len(osn_kontrol)
    
    if max_errors_mode == "гарантированно":
        max_errors_to_correct = k - 1
        print(f"Режим: гарантированное исправление до {max_errors_to_correct} ошибок")
    else:
        max_errors_to_correct = k
        print(f"Режим: попытка исправления до {max_errors_to_correct} ошибок (не гарантировано)")
    
    print(f"Количество оснований: {osn_kolvo}")
    print(f"Рабочий диапазон: {work_diapazon}")
    print(f"Количество контрольных оснований: {k}")
    print(f"Максимальное количество исправляемых ошибок: {max_errors_to_correct}")
    
    osn_int = [int(x) for x in osn]
    full_diapazon = prod(osn_int)
    chisl_int = [int(x) for x in chisl]
    
    for num_errors in range(1, max_errors_to_correct + 1):
        print(f"\nПроверка комбинаций из {num_errors} ошибок:")
        error_combinations = list(combinations(range(osn_kolvo), num_errors))
        combinations_count = len(error_combinations)
        print(f"Всего комбинаций: {combinations_count}")
        
        for error_positions in error_combinations:
            print(f"\nПроверка комбинации: позиции {error_positions}")
            print(f"Основания с ошибками: {[osn_int[pos] for pos in error_positions]}")
            
            # Создаем проекцию
            proj_osn = [osn_int[i] for i in range(osn_kolvo) if i not in error_positions]
            proj_chisl = [chisl_int[i] for i in range(osn_kolvo) if i not in error_positions]
            print(f"Проекционные основания: {proj_osn}")
            print(f"Проекционное число: {proj_chisl}")
            
            # Вычисляем ортогональные базисы для проекции
            B_proj = calculate_orthogonal_bases(proj_osn)
            if B_proj is None:
                continue
            
            proj_full_diapazon = prod(proj_osn)
            
            # Безопасное вычисление A_proj с модульной арифметикой
            A_proj = safe_dot_product(proj_chisl, B_proj, proj_full_diapazon)
            print(f"A_proj = {A_proj}")
            print(f"Сравнение: {A_proj} < {work_diapazon} = {A_proj < work_diapazon}")
            
            if A_proj < work_diapazon:
                print(f"Найдена корректная проекция! A_proj = {A_proj}")
                
                # Вычисляем полные базисы
                B_full = calculate_orthogonal_bases(osn_int)
                if B_full is None:
                    continue
                
                print(f"Полные ортогональные базисы: {B_full}")
                print('Создаем матрицу и вектор для решения системы уравнений')
                
                # Создаем матрицу системы
                M = []
                error_mods = [osn_int[pos] for pos in error_positions]
                print(f'Модули для ошибок: {error_mods}')
                
                for i, pos_i in enumerate(error_positions):
                    row = []
                    for j, pos_j in enumerate(error_positions):
                        value = B_full[pos_i] % error_mods[j]
                        row.append(value)
                        print(f'M[{i}][{j}] = {B_full[pos_i]} % {error_mods[j]} = {value}')
                    M.append(row)
                
                print(f"Матрица системы M:")
                for row in M:
                    print(row)
                
                # Создаем вектор b
                b = []
                
                # Вычисляем сумму для правильных позиций с модульной арифметикой
                s = 0
                for p in range(osn_kolvo):
                    if p not in error_positions:
                        term = (chisl_int[p] * B_full[p]) % full_diapazon
                        s = (s + term) % full_diapazon
                
                print(f"Сумма для правильных позиций s = {s}")
                
                for i, pos in enumerate(error_positions):
                    b_val = (A_proj - s) % osn_int[pos]
                    b.append(b_val)
                    print(f'b[{i}] = ({A_proj} - {s}) % {osn_int[pos]} = {b_val}')
                
                print(f"Вектор b: {b}")
                
                try:
                    print("Решение системы уравнений...")
                    solutions = solve_modular_system(M, b, error_mods)
                    
                    if solutions is None:
                        raise ValueError("Не удалось решить систему уравнений")
                    
                    print(f"Решения: {solutions}")
                    
                    # Применяем решения
                    print('Применяем решения')
                    corrected = chisl_int[:]
                    print('Исходное число: ', corrected)
                    
                    for i, pos in enumerate(error_positions):
                        corrected[pos] = solutions[i] % osn_int[pos]
                        print(f'corrected[{pos}] = {solutions[i]} % {osn_int[pos]} = {corrected[pos]}')
                    
                    print(f"Исправленное число: {corrected}")
                    
                    # Проверяем исправленное число с безопасным вычислением
                    A_corr = safe_dot_product(corrected, B_full, full_diapazon)
                    print(f"A_corr = {A_corr}")
                    print(f"Проверка исправления: {A_corr} < {work_diapazon} = {A_corr < work_diapazon}")
                    
                    if A_corr < work_diapazon:
                        corrections = [(pos, int(solutions[i])) for i, pos in enumerate(error_positions)]
                        print(f"Успешно исправлено {num_errors} ошибок!")
                        return corrections, f"Исправлено {num_errors} ошибок"
                    else:
                        print("Исправление не прошло проверку")
                        
                except Exception as e:
                    print(f"Ошибка при решении системы: {e}")
                    continue
    
    return None, "Не удалось обнаружить ошибки (возможно, их больше чем k)"

def main():    
    #Замер времени и памяти
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
    osn_inform = [int(input(f'Введите {i+1} информационное основание: ')) 
                  for i in range(osn_inform_kolvo)]
    osn_kontrol = [int(input(f'Введите {i+1} контрольное основание: ')) 
                   for i in range(osn_kontrol_kolvo)]
    osn = osn_inform + osn_kontrol
    print('\nВведите число для проверки:')
    chisl = [int(input(f'Введите остаток по основанию {osn[i]}: ')) 
             for i in range(len(osn))]
    
    print(f"Информационные основания: {osn_inform}")
    print(f"Контрольные основания: {osn_kontrol}")
    print(f"Все основания: {osn}")
    print(f"Проверяемое число: {chisl}")
    print(f"Режим работы: {max_errors_mode}")
    
    #Проверка числа методом ортогональных базисов
    print('Проверка методом ортогональных базисов')
    
    B = calculate_orthogonal_bases(osn)
    if B is None:
        print("Ошибка при вычислении ортогональных базисов")
        return
    
    full_diapazon = prod(osn)
    
    # Безопасное вычисление A с модульной арифметикой
    A = safe_dot_product(chisl, B, full_diapazon)
    
    work_diapazon = prod(osn_inform)
    print('Результат проверки')
    print(f"Вычисленное число A = {A}")
    print(f"Рабочий диапазон = {work_diapazon}")
    print(f"Проверка: {A} < {work_diapazon} = {A < work_diapazon}")
    
    if A < work_diapazon:
        print("Число корректно (в рабочем диапазоне)")
        
        #Замер ресурсов
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        execution_time = time.time() - start_time
        
        print(f"\nРесурсы программы:")
        print(f"Пиковое использование памяти: {peak / 1024:.2f} KB")
        print(f"Время выполнения: {execution_time:.6f} секунд")
        return
    
    print("Обнаружена ошибка (число вне рабочего диапазона)")
    
    #Применяем метод проекций
    corrections, message = check_with_projections(osn_inform, osn_kontrol, chisl, osn, max_errors_mode)
    
    #Вывод результатов
    print(f"Результат: {message}")
    
    if corrections:
        print("Рекомендуемые исправления:")
        corrected = chisl[:]
        for pos, correct_value in corrections:
            print(f"  - По основанию {osn[pos]}: {chisl[pos]} → {correct_value}")
            corrected[pos] = correct_value
        
        print(f"\nИсходное число: {chisl}")
        print(f"Исправленное число: {corrected}")
        
        # Проверка исправленного числа
        B_full = calculate_orthogonal_bases(osn)
        if B_full is not None:
            full_diapazon = prod(osn)
            A_corr = safe_dot_product(corrected, B_full, full_diapazon)
            print(f"\nПроверка исправленного числа: A = {A_corr}")
            print(f"Условие: {A_corr} < {work_diapazon} = {A_corr < work_diapazon}")
            
            if A_corr < work_diapazon:
                print("Исправление подтверждено!")
            else:
                print("Ошибка при исправлении!")
    
    #Замер ресурсов
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    execution_time = time.time() - start_time
    print(f"\nРесурсы программы:")
    print(f"Пиковое использование памяти: {peak / 1024:.2f} KB")
    print(f"Время выполнения: {execution_time:.6f} секунд")

if __name__ == "__main__":
    main()
