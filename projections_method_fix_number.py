import numpy
from itertools import combinations
from math import ceil
import time
import tracemalloc

#Вычисление ортогональных базисов
def calculate_orthogonal_bases(osn):
    print('Вычисление ортогональных базисов')
    full_diapazon = numpy.prod(osn.astype(object))
    print(f"1. Полный диапазон P = {full_diapazon}")

    P = full_diapazon // osn
    print(f"2. Величины P/pi:")
    for i in range(len(osn)):
        print(f"P[{i}] = {full_diapazon} // {osn[i]} = {P[i]}")
    beta = numpy.mod(P, osn)
    print(f"3. Величины βi = Pi mod pi:")
    for i in range(len(osn)):
        print(f"β[{i}] = {P[i]} mod {osn[i]} = {beta[i]}")
    print(f"4. Обратные элементы mi = inv(βi) mod pi:")
    m = numpy.zeros(len(osn), dtype=object)
    for i in range(len(osn)):
        try:
            m[i] = pow(int(beta[i]), -1, int(osn[i]))
            print(f"m[{i}] = inv({beta[i]}) mod {osn[i]} = {m[i]}")
        except ValueError:
            print(f"Ошибка: невозможно найти обратный элемент для β[{i}] = {beta[i]} mod {osn[i]}")
            return None
    result = m * P
    print(f"5. Ортогональные базисы Bi = mi * Pi:")
    for i in range(len(osn)):
        print(f"B[{i}] = {m[i]} * {P[i]} = {result[i]}")
    return result

#Вычисление методом проекций
def check_with_projections(osn_inform, osn_kontrol, chisl, osn, max_errors_mode):
    print('Вычисление методом проекций')
    osn_kolvo = len(osn)
    work_diapazon = numpy.prod(osn_inform.astype(object))
    k = len(osn_kontrol)
    #Определение максимального количества ошибок для исправления
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
    #Преобразуем osn в массив целых чисел
    osn_int = osn.astype(int)
    #Попытка исправления ошибок
    #Проверяются комбинации в зависимости от максимального количества исправляемых ошибок
    #Сначала все одинарные, затем двойные и т.д. До тех пор, пока не найдется ошибка (ошибки)
    for num_errors in range(1, max_errors_to_correct + 1):
        print(f"\nПроверка комбинаций из {num_errors} ошибок:")
        combinations_count = len(list(combinations(range(osn_kolvo), num_errors)))
        print(f"Всего комбинаций: {combinations_count}")
        for error_positions in combinations(range(osn_kolvo), num_errors):
            print(f"\nПроверка комбинации: позиции {error_positions}")
            print(f"Основания с ошибками: {osn[list(error_positions)]}")
            #Создаем проекцию
            proj_osn = numpy.delete(osn_int, error_positions)
            proj_chisl = numpy.delete(chisl, error_positions)
            print(f"Проекционные основания: {proj_osn}")
            print(f"Проекционное число: {proj_chisl}")            
            # Вычисляем ортогональные базисы для проекции
            B_proj = calculate_orthogonal_bases(proj_osn)
            if B_proj is None:
                continue
            A_proj = numpy.mod(numpy.sum(proj_chisl * B_proj), numpy.prod(proj_osn.astype(object)))
            print(f"A_proj = {A_proj}")
            print(f"Сравнение: {A_proj} < {work_diapazon} = {A_proj < work_diapazon}")
            if A_proj < work_diapazon:
                print(f"Найдена корректная проекция! A_proj = {A_proj}")
                # Вычисляем полные базисы
                B_full = calculate_orthogonal_bases(osn_int)
                if B_full is None:
                    continue
                print(f"Полные ортогональные базисы: {B_full}")
                
                #Создаем матрицу и вектор для решения системы уравнений
                print('Создаем матрицу и вектор для решения системы уравнений')
                M = numpy.zeros((num_errors, num_errors), dtype=object)
                print('M = ', M)
                for i, pos in enumerate(error_positions):
                    for j, p in enumerate(error_positions):
                        M[i,j] = B_full[pos] % osn_int[p]
                        print('M[', i, '][', j, '] = ', B_full[pos], ' % ', osn_int[p])
                        print('M[', i, '][', j, '] = ', M[i, j])
                print(f"Матрица системы M:")
                for i in range(num_errors):
                    print(f"{M[i]}")
                b = numpy.zeros(num_errors, dtype=object)
                print('b = ', b)
                for i, pos in enumerate(error_positions):
                    # Вычисляем сумму для правильных позиций
                    s = 0
                    for p in range(osn_kolvo):
                        if p not in error_positions:
                            s += chisl[p] * B_full[p]
                    s %= numpy.prod(osn_int.astype(object))
                    print(s)
                    
                    b[i] = (A_proj - s) % osn_int[pos]
                print(f"Вектор b: {b}")
                try:
                    print("Решение системы уравнений...")
                    #Используем точное решение для больших чисел
                    solutions = numpy.zeros(num_errors, dtype=object)
                    print('solutions = ', solutions)
                    #Для одиночной ошибки - простое решение
                    if num_errors == 1:
                        solutions[0] = (b[0] * pow(int(M[0,0]), -1, int(osn_int[error_positions[0]]))) % osn_int[error_positions[0]]
                        print(solutions[0])
                    else:
                        #Для множественных ошибок используем метод Крамера или Гаусса
                        #Упрощенная реализация для демонстрации
                        det_M = numpy.linalg.det(M.astype(float))
                        print('det_M = ', det_M)
                        if abs(det_M) < 1e-10:
                            raise numpy.linalg.LinAlgError("Матрица вырождена")
                        for i in range(num_errors):
                            M_temp = M.copy().astype(float)
                            M_temp[:, i] = b.astype(float)
                            det_temp = numpy.linalg.det(M_temp)
                            solutions[i] = round(det_temp / det_M) % osn_int[error_positions[i]]
                    
                    print(f"Решения: {solutions}")
                    
                    #Применяем решения
                    print('Применяем решения')
                    corrected = chisl.copy()
                    print('corrected = ', corrected)
                    for i, pos in enumerate(error_positions):
                        corrected[pos] = solutions[i] % osn_int[pos]
                        print('corrected[', pos, '] = ', corrected[pos])
                    
                    print(f"Исправленное число: {corrected}")
                    
                    #Проверяем исправленное число
                    A_corr = numpy.mod(numpy.sum(corrected * B_full), numpy.prod(osn_int.astype(object)))
                    print(f"A_corr = {A_corr}")
                    print(f"Проверка исправления: {A_corr} < {work_diapazon} = {A_corr < work_diapazon}")
                    
                    if A_corr < work_diapazon:
                        corrections = [(pos, int(solutions[i])) for i, pos in enumerate(error_positions)]
                        print(f"Успешно исправлено {num_errors} ошибок!")
                        return corrections, f"Исправлено {num_errors} ошибок"
                    else:
                        print("Исправление не прошло проверку")
                        
                except (numpy.linalg.LinAlgError, ValueError) as e:
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
    osn_inform = numpy.array([int(input(f'Введите {i+1} информационное основание: ')) 
                            for i in range(osn_inform_kolvo)], dtype=object)
    osn_kontrol = numpy.array([int(input(f'Введите {i+1} контрольное основание: ')) 
                            for i in range(osn_kontrol_kolvo)], dtype=object)
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    print('\nВведите число для проверки:')
    chisl = numpy.array([int(input(f'Введите остаток по основанию {osn[i]}: ')) 
                       for i in range(len(osn))], dtype=object)
    
    print(f"Информационные основания: {osn_inform}")
    print(f"Контрольные основания: {osn_kontrol}")
    print(f"Все основания: {osn}")
    print(f"Проверяемое число: {chisl}")
    print(f"Режим работы: {max_errors_mode}")
    
    #Проверка числа методом ортогональных базисов
    print('Проверка методом ортогональных базисов')
    
    B = calculate_orthogonal_bases(osn.astype(int))
    if B is None:
        print("Ошибка при вычислении ортогональных базисов")
        return
    
    A = numpy.mod(numpy.sum(chisl * B), numpy.prod(osn.astype(object)))
    work_diapazon = numpy.prod(osn_inform.astype(object))
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
        corrected = chisl.copy()
        for pos, correct_value in corrections:
            print(f"  - По основанию {osn[pos]}: {chisl[pos]} → {correct_value}")
            corrected[pos] = correct_value
        
        print(f"\nИсходное число: {chisl}")
        print(f"Исправленное число: {corrected}")
        
        # Проверка исправленного числа
        B_full = calculate_orthogonal_bases(osn.astype(int))
        if B_full is not None:
            A_corr = numpy.mod(numpy.sum(corrected * B_full), numpy.prod(osn.astype(object)))
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

main()
