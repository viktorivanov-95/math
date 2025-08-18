import numpy
import time
import tracemalloc

def check_number_with_orthogonal_bases(osn_inform, osn_kontrol, chisl):
    #Подсчет затрат памяти и времени
    start_time = time.time()
    tracemalloc.start()
    
    print("\n1. Объединение информационных и контрольных оснований:")
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    osn_kolvo = len(osn)
    print(f"Основания системы: {osn}")
    print(f"Введенное число: {chisl}")

    #2. Вычисление ортогональных базисов
    print("\n2. Вычисление ортогональных базисов:")
    #2.1 Поиск рабочего и полного диапазонов
    print("\n2.1 Поиск рабочего R и полного P диапазонов:")
    work_diapazon = numpy.prod(osn_inform.astype(object))
    full_diapazon = numpy.prod(osn.astype(object))
    print(f"Рабочий диапазон R (произведение информационных оснований) = {work_diapazon}")
    print(f"Полный диапазон P (произведение всех оснований) = {full_diapazon}")
    #2.2 Поиск величины P/pi = Pi для каждого основания
    print("\n2.2 Поиск величины P/pi = Pi для каждого основания:")
    P = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        P[i] = full_diapazon // osn[i]
        print(f"P[{i}] = {full_diapazon} // {osn[i]} = {P[i]}")
    print(f"Величина P = {P}")
    #2.3 Поиск величины βi = Pi(modpi)
    print("\n2.3 Поиск величины βi = Pi(modpi):")
    beta = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        beta[i] = P[i] % osn[i]
        print(f"β[{i}] = {P[i]} mod {osn[i]} = {beta[i]}")
    print(f"Величина β = {beta}")
    #2.4 Нахождение веса mi базисов
    print("\n2.4 Нахождение весов mi базисов (обратный элемент βi по модулю pi):")
    m = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        try:
            m[i] = pow(int(beta[i]), -1, int(osn[i]))
            print(f"m[{i}] = inv({beta[i]}) mod {osn[i]} = {m[i]}")
        except ValueError:
            print(f"Ошибка: невозможно найти обратный элемент для beta[{i}] = {beta[i]} по модулю {osn[i]}")
            return None, "Ошибка при вычислении обратного элемента"
    print(f"Веса базисов m = {m}")
    #2.5 Вычисление ортогональных базисов системы
    print("\n2.5 Вычисление ортогональных базисов системы Bi = mi * Pi:")
    B = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        B[i] = m[i] * P[i]
        print(f"B[{i}] = {m[i]} * {P[i]} = {B[i]}")
    print(f"Ортогональные базисы системы B = {B}")
    #2.6 Нахождение числа A в десятичной системе
    print("\n2.6 Нахождение числа A в десятичной системе:")
    A = 0
    for i in range(osn_kolvo):
        old_A = A
        A += chisl[i] * B[i]
        print(f"A[{i}] = {old_A} + {chisl[i]} * {B[i]} = {A}")

    new_A = A
    A = new_A % full_diapazon
    print(f"\nA = {new_A} mod {full_diapazon} = {A}")
    #2.7 Проверка числа
    print("\n2.7 Проверка числа:")
    print(f"Сравниваем A = {A} с рабочим диапазоном R = {work_diapazon}:")
    if A > work_diapazon:
        result = f"{chisl} - представление числа, содержащее ошибку (ошибки)"
    else:
        result = f"{chisl} - число не содержит ошибки (ошибок)"
    print(result)
    
    #Замер памяти и времени выполнения
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    execution_time = time.time() - start_time
    print("\nМетрики выполнения:")
    print(f"Пиковое использование памяти: {peak / 1024:.2f} KB")
    print(f"Время выполнения: {execution_time:.6f} секунд")
    
    return A, result

#Ввод данных
print("Основания нужно вводить в порядке возрастания: сначала информационные, затем контрольные.")
print("Они должны быть взаимно простыми.")
osn_inform_kolvo = int(input('Введите количество информационных оснований: '))
osn_kontrol_kolvo = int(input('Введите количество контрольных оснований: '))
osn_inform = numpy.zeros((osn_inform_kolvo), dtype=object)
osn_kontrol = numpy.zeros((osn_kontrol_kolvo), dtype=object)
for i in range(osn_inform_kolvo):
    osn_inform[i] = int(input(f'Введите {i+1} информационное основание: '))
for i in range(osn_kontrol_kolvo):
    osn_kontrol[i] = int(input(f'Введите {i+1} контрольное основание: '))
#Ввод числа для проверки
osn_kolvo = osn_inform_kolvo + osn_kontrol_kolvo
chisl = numpy.zeros((osn_kolvo), dtype=object)
for i in range(osn_kolvo):
    chisl[i] = int(input(f'Введите {i+1} остаток числа: '))
#Проверка числа
A, result = check_number_with_orthogonal_bases(osn_inform, osn_kontrol, chisl)
if A is not None:
    print("\nИтоговый результат:")
    print(f"Вычисленное десятичное число A = {A}")
    print(f"Результат проверки: {result}")
