import numpy
import time
import tracemalloc
import math
from decimal import Decimal, getcontext, ROUND_DOWN, ROUND_HALF_UP

#Метод ортогональных базисов
def check_with_orthogonal_bases(osn_inform, osn_kontrol, chisl):
    print("\nМетод ортогональных базисов")
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    osn_kolvo = len(osn)
    work_diapazon = numpy.prod(osn_inform.astype(object))
    full_diapazon = numpy.prod(osn.astype(object))
    print(f"Основания системы: {osn}")
    print(f"Остатки числа: {chisl}")
    print(f"Рабочий диапазон (R): {work_diapazon}")
    print(f"Полный диапазон (P): {full_diapazon}")
    # 1. Вычисление Pi = P/pi
    print("\n1. Вычисление Pi = P/pi для каждого основания")
    P = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        P[i] = full_diapazon // osn[i]
        print(f"P[{i}] = {full_diapazon} // {osn[i]} = {P[i]}")
    # 2. Вычисление βi = Pi mod pi
    print("\n2. Вычисление βi = Pi mod pi")
    beta = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        beta[i] = P[i] % osn[i]
        print(f"β[{i}] = {P[i]} mod {osn[i]} = {beta[i]}")
    # 3. Вычисление mi (обратные к βi по модулю pi)
    print("\n3. Вычисление весов mi (обратные к βi)")
    m = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        m[i] = pow(int(beta[i]), -1, int(osn[i]))
        print(f"m[{i}] = inv({beta[i]}) mod {osn[i]} = {m[i]}")
    # 4. Вычисление ортогональных базисов Bi = mi * Pi
    print("\n4. Вычисление ортогональных базисов Bi")
    B = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        B[i] = m[i] * P[i]
        print(f"B[{i}] = {m[i]} * {P[i]} = {B[i]}")
    # 5. Вычисление числа A
    print("\n5. Вычисление числа A")
    A = 0
    for i in range(osn_kolvo):
        if i == 0:
            print(f"A[{i}] = 0 + {int(chisl[i])} * {int(B[i])} = {int(chisl[i]*B[i])}")
        else:
            print(f"A[{i}] = {int(A)} + {int(chisl[i])} * {int(B[i])} = {int(A + chisl[i]*B[i])}")
        A = int(A + chisl[i] * B[i])
    A = numpy.mod(A, full_diapazon)
    print(f"\nA = {A} mod {full_diapazon} = {A}")
    # Проверка на ошибки
    print("\nПроверка результата:")
    if A > work_diapazon:
        print(f"Ошибка: число {A} находится вне рабочего диапазона (0-{work_diapazon-1})")
    else:
        print(f"Число {A} находится в рабочем диапазоне (0-{work_diapazon-1})")
    
    return A

#Метод перевода числа в ОПСС
def check_with_opss(osn_inform, osn_kontrol, chisl):
    print("\nМетод перевода числа в ОПСС")
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    osn_kolvo = len(osn)
    print(f"Основания системы: {osn}")
    print(f"Остатки числа: {chisl}")
    #1. Вычисление констант c_ij (обратные элементы)
    print("\n1. Вычисление констант c_ij (обратные элементы)")
    c = numpy.zeros((osn_kolvo, osn_kolvo), dtype=object)
    for i in range(osn_kolvo):
        for j in range(i + 1, osn_kolvo):
            c[i][j] = pow(int(osn[i]), -1, int(osn[j]))
            print(f"c[{i}][{j}] = inv({osn[i]}) mod {osn[j]} = {c[i][j]}")
    #2. Перевод в ОПСС
    print("\n2. Перевод числа в ОПСС")
    opss = numpy.zeros(osn_kolvo, dtype=object)
    opss[0] = chisl[0] % osn[0]
    print(f"opss[0] = {chisl[0]} mod {osn[0]} = {opss[0]}")
    #Последующие разряды
    for i in range(1, osn_kolvo):
        temp = chisl[i]
        for j in range(i):
            temp = (temp - opss[j]) * c[j][i] % osn[i]
        opss[i] = temp % osn[i]
        print(f"opss[{i}] = {temp} mod {osn[i]} = {opss[i]}")
    print(f"\nЧисло в ОПСС: {opss}")
    #3. Проверка контрольных остатков
    print("\n3. Проверка контрольных остатков:")
    error = False
    for i in range(len(osn_inform), osn_kolvo):
        if opss[i] != 0:
            print(f"Ошибка: opss[{i}] = {opss[i]} ≠ 0 (по контрольному основанию {osn[i]})")
            error = True
    if not error:
        print("Все контрольные остатки равны 0 - ошибок не обнаружено")
    #4. Перевод в десятичную систему
    print("\n4. Перевод из ОПСС в десятичную систему")
    A_opss = 0
    product = 1
    for i in range(osn_kolvo):
        term = opss[i] * product
        print(f"{opss[i]} * {product} = {term}")
        A_opss += term
        if i < osn_kolvo - 1:
            product *= osn[i]
    print(f"\nИтоговое число: {A_opss}")
    return A_opss

#Проверка методом совместного использования ортогональных базисов и ОПСС
def check_with_combined_method(osn_inform, osn_kontrol, chisl):
    print("\nКомбинированный метод")
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    osn_kolvo = len(osn)
    full_diapazon = numpy.prod(osn.astype(object))
    print(f"Основания системы: {osn}")
    print(f"Остатки числа: {chisl}")
    print(f"Полный диапазон (P): {full_diapazon}")
    #1. Вычисление ортогональных базисов
    print("\n1. Вычисление ортогональных базисов")
    P = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        P[i] = full_diapazon // osn[i]
        print(f"P[{i}] = {full_diapazon} // {osn[i]} = {P[i]}")
    beta = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        beta[i] = P[i] % osn[i]
        print(f"β[{i}] = {P[i]} mod {osn[i]} = {beta[i]}")
    m = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        m[i] = pow(int(beta[i]), -1, int(osn[i]))
        print(f"m[{i}] = inv({beta[i]}) mod {osn[i]} = {m[i]}")
    B = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        B[i] = m[i] * P[i]
        print(f"B[{i}] = {m[i]} * {P[i]} = {B[i]}")
    #2. Перевод базисов в СОК
    print("\n2. Перевод базисов в СОК")
    B_modular = numpy.zeros((osn_kolvo, osn_kolvo), dtype=object)
    for i in range(osn_kolvo):
        B_modular[i][:] = [B[j] % osn[i] for j in range(osn_kolvo)]
        print(f"B_modular[{i}] = B mod {osn[i]} = {B_modular[i][:]}")
    #3. Вычисление констант для ОПСС
    print("\n3. Вычисление констант c_ij для ОПСС")
    c = numpy.zeros((osn_kolvo, osn_kolvo), dtype=object)
    for i in range(osn_kolvo):
        for j in range(i + 1, osn_kolvo):
            c[i][j] = pow(int(osn[i]), -1, int(osn[j]))
            print(f"c[{i}][{j}] = inv({osn[i]}) mod {osn[j]} = {c[i][j]}")
    #4. Разложение базисов в ОПСС
    print("\n4. Разложение базисов в ОПСС")
    B_opss = numpy.zeros((osn_kolvo, osn_kolvo), dtype=object)
    for k in range(osn_kolvo):
        print(f"\nРазложение базиса B[{k}] = {B[k]}:")
        B_opss[k][0] = B_modular[k][0] % osn[0]
        print(f"B_opss[{k}][0] = {B_modular[k][0]} mod {osn[0]} = {B_opss[k][0]}")
        for i in range(1, osn_kolvo):
            temp = B_modular[k][i]
            for j in range(i):
                temp = (temp - B_opss[k][j]) * c[j][i] % osn[i]
            B_opss[k][i] = temp % osn[i]
            print(f"B_opss[{k}][{i}] = {temp} mod {osn[i]} = {B_opss[k][i]}")
    #5. Расчет числа в ОПСС
    print("\n5. Расчет числа в ОПСС")
    B_opss_chisl = numpy.zeros((osn_kolvo, osn_kolvo), dtype=object)
    for i in range(osn_kolvo):
        for j in range(osn_kolvo):
            B_opss_chisl[i][j] = chisl[i] * B_opss[i][j]
            print(f"B_opss_chisl[{i}][{j}] = {chisl[i]} * {B_opss[i][j]} = {B_opss_chisl[i][j]}")
    opss = numpy.zeros(osn_kolvo, dtype=object)
    for j in range(osn_kolvo):
        for i in range(osn_kolvo):
            opss[j] += B_opss_chisl[i][j]
        print(f"opss[{j}] = сумма по столбцу {j} = {opss[j]}")
    #6. Учет переносов
    print("\n6. Учет переносов")
    opss_total = numpy.zeros(osn_kolvo, dtype=object)
    perenos = 0
    for i in range(osn_kolvo):
        total = opss[i] + perenos
        opss_total[i] = total % osn[i]
        new_perenos = total // osn[i]
        print(f"opss_total[{i}] = ({opss[i]} + {perenos}) mod {osn[i]} = {opss_total[i]}")
        print(f"perenos = ({opss[i]} + {perenos}) // {osn[i]} = {new_perenos}")
        perenos = new_perenos
    print(f"\nЧисло в ОПСС с учетом переносов: {opss_total}")
    #7. Перевод в десятичную систему
    print("\n7. Перевод в десятичную систему")
    A_combined = 0
    product = 1
    for i in range(osn_kolvo):
        term = opss_total[i] * product
        print(f"{opss_total[i]} * {product} = {term}")
        A_combined += term
        if i < osn_kolvo - 1:
            product *= osn[i]
    print(f"\nИтоговое число: {A_combined}")
    return A_combined

#Метод относительных величин с поддержкой больших чисел
def relative_values_method(osn_inform, osn_kontrol, chisl, precision=10):
    print(f"\nМетод относительных величин (точность: {precision} знаков)")
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    osn_kolvo = len(osn)
    #Используем Decimal для вычисления диапазонов
    work_diapazon = Decimal(1)
    for base in osn_inform:
        work_diapazon *= Decimal(str(base))
    full_diapazon = work_diapazon
    for base in osn_kontrol:
        full_diapazon *= Decimal(str(base))
    print(f"Основания системы: {osn}")
    print(f"Остатки числа: {chisl}")
    print(f"Рабочий диапазон (R): {work_diapazon}")
    print(f"Полный диапазон (P): {full_diapazon}")
    #Настройка точности Decimal
    getcontext().prec = precision + 50
    getcontext().rounding = ROUND_HALF_UP
    #1. Вычисление Pi = P/pi
    print("\n1. Вычисление Pi = P/pi для каждого основания")
    P = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        base_decimal = Decimal(str(osn[i]))
        P_decimal = full_diapazon / base_decimal
        # Ограничиваем количество знаков после запятой
        P_decimal = P_decimal.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
        P[i] = P_decimal
        print(f"P[{i}] = {full_diapazon} / {base_decimal} = {P[i]}")
    #2. Вычисление коэффициентов k[i]
    print("\n2. Вычисление коэффициентов k[i] = (Pi^(-1) mod pi) / pi")
    k = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        # Вычисляем обратный элемент с целыми числами
        Pi_int = int(Decimal(P[i]).to_integral_value(rounding=ROUND_HALF_UP))
        base_int = int(osn[i])
        inv_Pi = pow(Pi_int, -1, base_int)
        # Вычисляем k[i] с ограниченной точностью
        k_decimal = Decimal(str(inv_Pi)) / Decimal(str(base_int))
        k_decimal = k_decimal.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
        k[i] = k_decimal
        print(f"k[{i}] = inv({Pi_int}) mod {base_int} / {base_int} = {k[i]}")
    #3. Вычисление A/P
    print("\n3. Вычисление A/P = Σ(chisl[i] * k[i])")
    A_relative = Decimal('0')
    for i in range(osn_kolvo):
        term = Decimal(str(chisl[i])) * k[i]
        term = term.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
        if i == 0:
            print(f"A/P[{i}] = 0 + {chisl[i]} * {k[i]} = {term}")
        else:
            print(f"A/P[{i}] = {A_relative} + {chisl[i]} * {k[i]} = {A_relative + term}")
        A_relative += term
        A_relative = A_relative.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
    #4. Вычисление A
    print(f"\n4. Вычисление A = (A/P mod 1) * P")
    fractional_part = A_relative % Decimal('1')
    fractional_part = fractional_part.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
    print(f"Дробная часть A/P = {fractional_part}")
    A_relative = fractional_part * full_diapazon
    A_relative = A_relative.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
    print(f"A = {fractional_part} * {full_diapazon} = {A_relative}")
    #5. Проверка на ошибки
    print("\n5. Проверка результата:")
    if A_relative > work_diapazon:
        print(f"Ошибка: число {A_relative} находится вне рабочего диапазона (0-{work_diapazon-1})")
    else:
        print(f"Число {A_relative} находится в рабочем диапазоне (0-{work_diapazon-1})")
    
    return A_relative

#Начало выполнения программы
start_time = time.time()
tracemalloc.start()

print("Основания нужно вводить в порядке возрастания: сначала информационные, затем контрольные.")
print("Они должны быть взаимно простыми.")

#Ввод данных
osn_inform_kolvo = int(input('Введите количество информационных оснований: '))
osn_kontrol_kolvo = int(input('Введите количество контрольных оснований: '))
osn_inform = numpy.zeros(osn_inform_kolvo, dtype=object)
osn_kontrol = numpy.zeros(osn_kontrol_kolvo, dtype=object)

for i in range(osn_inform_kolvo):
    osn_inform[i] = int(input(f'Введите {i+1} информационное основание: '))
for i in range(osn_kontrol_kolvo):
    osn_kontrol[i] = int(input(f'Введите {i+1} контрольное основание: '))

osn = numpy.concatenate((osn_inform, osn_kontrol))
osn_kolvo = len(osn)
chisl = numpy.zeros(osn_kolvo, dtype=object)

for i in range(osn_kolvo):
    chisl[i] = int(input(f'Введите остаток по основанию {osn[i]}: '))

#Настройка точности метода относительных величин
print("\nНастройка точности метода относительных величин")
while True:
    try:
        precision = int(input('Введите количество знаков после запятой для точности (1-50): '))
        if 1 <= precision <= 50:
            break
        else:
            print("Ошибка: введите число от 1 до 50")
    except ValueError:
        print("Ошибка: введите целое число")

#Вывод введенных данных
print("\nВведенные данные")
print('Информационные основания:', osn_inform)
print('Контрольные основания:', osn_kontrol)
print('Остатки числа:', chisl)

#Вычисление диапазонов с использованием Decimal для больших чисел
work_diapazon_dec = Decimal(1)
for base in osn_inform:
    work_diapazon_dec *= Decimal(str(base))

full_diapazon_dec = work_diapazon_dec
for base in osn_kontrol:
    full_diapazon_dec *= Decimal(str(base))

print('Рабочий диапазон R:', work_diapazon_dec)
print('Полный диапазон P:', full_diapazon_dec)

print("Проверка числа всеми методами")

#Метод относительных величин
A_relative = relative_values_method(osn_inform, osn_kontrol, chisl, precision)
# Точные методы
A_orth = check_with_orthogonal_bases(osn_inform, osn_kontrol, chisl)
A_opss = check_with_opss(osn_inform, osn_kontrol, chisl)
A_combined = check_with_combined_method(osn_inform, osn_kontrol, chisl)
#Сравнение результатов
print("Сравнение результатов")
print(f"Метод относительных величин: A = {A_relative}")
print(f"Метод ортогональных базисов: A = {A_orth}")
print(f"Метод ОПСС: A = {A_opss}")
print(f"Комбинированный метод: A = {A_combined}")
#Преобразуем точные результаты в Decimal для сравнения
A_orth_dec = Decimal(str(A_orth))
A_opss_dec = Decimal(str(A_opss))
A_combined_dec = Decimal(str(A_combined))

print(f"\nРазница с методом ортогональных базисов: {abs(A_relative - A_orth_dec)}")
print(f"Разница с методом ОПСС: {abs(A_relative - A_opss_dec)}")
print(f"Разница с комбинированным методом: {abs(A_relative - A_combined_dec)}")

if A_orth == A_opss == A_combined:
    print("\nТочные методы дали одинаковый результат")
else:
    print("\nТочные методы дали разные результаты!")

#Проверка на ошибки
print("Проверка на ошибки")
print(f'Метод относительных величин: число {A_relative}', 
      'содержит ошибки' if A_relative > work_diapazon_dec else 'не содержит ошибок')

print(f'Точные методы: число {A_orth}', 
      'содержит ошибки' if A_orth > work_diapazon_dec else 'не содержит ошибок')

#Замер ресурсов
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
execution_time = time.time() - start_time

print(f'Пиковое использование памяти: {peak / 1024:.2f} KB')
print(f'Время выполнения: {execution_time:.6f} секунд')
print(f'Точность вычислений: {precision} знаков после запятой')
