import numpy
import time
import tracemalloc

def optimized_opss_conversion(osn_inform, osn_kontrol, chisl):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    osn_kolvo = len(osn)
    
    #Поиск констант c_ij:
    c = numpy.zeros((osn_kolvo-1, osn_kolvo), dtype=object)
    
    #Вычисляем только нужные константы (нижняя треугольная матрица)
    for j in range(osn_kolvo-1):
        for i in range(j+1):
            try:
                c[j][i] = pow(int(osn[i]), -1, int(osn[j+1]))
                print(f'c[{j}][{i}] = {c[j][i]} (поскольку {c[j][i]} * {osn[i]} ≡ 1 mod {osn[j+1]})')
            except ValueError:
                print(f'Невозможно найти c[{j}][{i}] - обратный элемент для {osn[i]} mod {osn[j+1]}')
                return None
    
    #Матрица констант c_ij:
    print(c[:,:osn_kolvo-1])
    
    #Расчет числа в ОПСС:
    opss = numpy.zeros(osn_kolvo, dtype=object)
    opss[0] = chisl[0] % osn[0]
    print(f'opss[0] = {chisl[0]} mod {osn[0]} = {opss[0]}')
    
    if osn_kolvo > 1:
        temp = (chisl[1] - opss[0]) * c[0][0]
        opss[1] = temp % osn[1]
        print(f'opss[1] = ({chisl[1]} - {opss[0]}) * {c[0][0]} = {temp} mod {osn[1]} = {opss[1]}')
    
    for i in range(2, osn_kolvo):
        #Первый шаг вычисления
        opss[i] = (chisl[i] - opss[0]) * c[i-1][0]
        print(f'\nНачальное opss[{i}] = ({chisl[i]} - {opss[0]}) * {c[i-1][0]} = {opss[i]}')
        
        #Коррекция с использованием предыдущих значений
        for j in range(1, i):
            opss[i] = (opss[i] - opss[j]) * c[i-1][j]
            print(f'Шаг {j}: opss[{i}] = (пред.результат - {opss[j]}) * {c[i-1][j]} = {opss[i]}')
        
        opss[i] %= osn[i]
        print(f'Финальное opss[{i}] = {opss[i]} mod {osn[i]}')
    
    print('\nЧисло в ОПСС равно:', opss)
    return opss

#Проверка методом ортогональных базисов
def check_with_orthogonal_bases(osn_inform, osn_kontrol, chisl):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    osn_kolvo = len(osn)
    full_diapazon = numpy.prod(osn.astype(object))
    
    print('\n1. Поиск величин P/pi = Pi для каждого основания:')
    P = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        P[i] = full_diapazon // osn[i]
        print(f'P[{i}] = {full_diapazon} / {osn[i]} = {P[i]}')
        
    print('\n2. Поиск величин βi = Pi mod pi:')
    beta = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        beta[i] = P[i] % osn[i]
        print(f'β[{i}] = {P[i]} mod {osn[i]} = {beta[i]}')
    
    print('\n3. Нахождение весов mi базисов (mi * βi ≡ 1 mod pi):')
    m = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        try:
            m[i] = pow(int(beta[i]), -1, int(osn[i]))
            print(f'm[{i}] = {beta[i]}^(-1) mod {osn[i]} = {m[i]} (поскольку {m[i]} * {beta[i]} ≡ 1 mod {osn[i]})')
        except ValueError:
            print(f'Yевозможно найти обратный элемент для β[{i}] = {beta[i]} по модулю {osn[i]}')
            return None
    
    print('\n4. Вычисление ортогональных базисов системы Bi = mi * Pi:')
    B = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        B[i] = m[i] * P[i]
        print(f'B[{i}] = {m[i]} * {P[i]} = {B[i]}')
    
    print('\n5. Нахождение числа A в позиционной системе счисления:')
    A = 0
    for i in range(osn_kolvo):
        prev_A = A
        A += chisl[i] * B[i]
        print(f'A[{i}] = {prev_A} + {chisl[i]} * {B[i]} = {A}')
    
    A = A % full_diapazon
    print(f'\nИтоговое A = {A} mod {full_diapazon} = {A}')
    return A

#Начало программы
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

print('\nВведите остатки числа по каждому основанию:')
chisl = numpy.zeros(len(osn_inform) + len(osn_kontrol), dtype=object)
for i in range(len(osn_inform)):
    chisl[i] = int(input(f'Введите остаток по информационному основанию {osn_inform[i]}: '))
for i in range(len(osn_kontrol)):
    chisl[len(osn_inform)+i] = int(input(f'Введите остаток по контрольному основанию {osn_kontrol[i]}: '))

osn = numpy.concatenate((osn_inform, osn_kontrol))
work_diapazon = numpy.prod(osn_inform.astype(object))
full_diapazon = numpy.prod(osn.astype(object))

print('\nВведенные данные:')
print('Информационные основания:', osn_inform)
print('Контрольные основания:', osn_kontrol)
print('Все основания:', osn)
print('Остатки числа:', chisl)
print(f'Рабочий диапазон R = {work_diapazon}')
print(f'Полный диапазон P = {full_diapazon}')

#Ппроверка методом ОПСС
print('Перевод числа из СОК в ОПСС')
opss = optimized_opss_conversion(osn_inform, osn_kontrol, chisl)

if opss is not None:
    #Проверка на ошибки по контрольным основаниям
    print('\nПроверка числа на наличие ошибок:')
    error_found = False
    for i in range(osn_inform_kolvo, len(osn)):
        if opss[i] != 0:
            print(f'Обнаружена ошибка: opss[{i}] = {opss[i]} ≠ 0 (по контрольному основанию {osn[i]})')
            error_found = True

    if not error_found:
        print('Ошибок не обнаружено: все контрольные остатки в ОПСС равны 0')

    #Перевод ОПСС в десятичное число
    print('\nПеревод числа из ОПСС в десятичную систему:')
    opss_check = 0
    for i in range(len(osn)):
        product = numpy.prod(osn[:i].astype(object)) if i > 0 else 1
        term = opss[i] * product
        print(f'{opss[i]} * {product} = {term}')
        opss_check += term
    print(f'Итоговое число: {opss_check}')

#Проверка методом ортогональных базисов
print('Проверка методом ортогональных базисов')
A_orth = check_with_orthogonal_bases(osn_inform, osn_kontrol, chisl)

#Сравнение результатов
if opss is not None and A_orth is not None:
    print('Сравнение результатов')
    print(f'Число из ОПСС: {opss_check}')
    print(f'Число из ортогональных базисов: {A_orth}')

    if opss_check == A_orth:
        print('Результаты совпадают')
    else:
        print('Результаты не совпадают!')
    
    #Проверка на ошибки по рабочему диапазону
    if A_orth > work_diapazon:
        print(f'Обнаружена ошибка: число {A_orth} находится вне рабочего диапазона (0-{work_diapazon-1})')
    else:
        print(f'Число {A_orth} находится в рабочем диапазоне (0-{work_diapazon-1})')

#Замер памяти и времени
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
execution_time = time.time() - start_time

print('Ресурсы программы:')
print(f'Использовано памяти: {current/1024:.2f} KB')
print(f'Пиковое использование памяти: {peak/1024:.2f} KB')
print(f'Время выполнения: {execution_time * 1e6:.2f} микросекунд')
