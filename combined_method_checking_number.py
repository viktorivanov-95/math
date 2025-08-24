import numpy
import time
import tracemalloc

def check_with_orthogonal_bases(osn_inform, osn_kontrol, chisl):
    print("\nПроверка методом ортогональных базисов")
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    osn_kolvo = len(osn)
    work_diapazon = numpy.prod(osn_inform)
    full_diapazon = numpy.prod(osn)
    
    print(f"Основания системы: {osn}")
    print(f"Остатки числа: {chisl}")
    print(f"Рабочий диапазон (R): {work_diapazon}")
    print(f"Полный диапазон (P): {full_diapazon}")

    #1. Вычисление Pi = P/pi
    print("Вычисление Pi = P/pi для каждого основания")
    P = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        P[i] = full_diapazon // osn[i]
        print(f"P[{i}] = {full_diapazon} // {osn[i]} = {P[i]}")
    #2. Вычисление βi = Pi mod pi
    print("Вычисление βi = Pi mod pi")
    beta = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        beta[i] = P[i] % osn[i]
        print(f"β[{i}] = {P[i]} mod {osn[i]} = {beta[i]}")
    #3. Вычисление mi (обратные к βi по модулю pi)
    print("Вычисление весов mi (обратные к βi)")
    m = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        try:
            m[i] = pow(int(beta[i]), -1, int(osn[i]))
            print(f"m[{i}] = inv({beta[i]}) mod {osn[i]} = {m[i]}")
        except ValueError:
            print(f"Ошибка: невозможно найти обратный элемент для β[{i}]")
            return None
    #4. Вычисление ортогональных базисов Bi = mi * Pi
    print("Вычисление ортогональных базисов Bi")
    B = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        B[i] = m[i] * P[i]
        print(f"B[{i}] = {m[i]} * {P[i]} = {B[i]}")
    #5. Вычисление числа A
    print("Вычисление числа A")
    A = 0
    for i in range(osn_kolvo):
        if i == 0:
            print(f"A[{i}] = 0 + {int(chisl[i])} * {int(B[i])} = {int(chisl[i]*B[i])}")
        else:
            print(f"A[{i}] = {int(A)} + {int(chisl[i])} * {int(B[i])} = {int(A + chisl[i]*B[i])}")
        A = int(A + chisl[i] * B[i])
    
    A = numpy.mod(A, full_diapazon)
    print(f"\nA = {A} mod {full_diapazon} = {A}")
    
    #Проверка на ошибки
    print("\nПроверка результата:")
    if A > work_diapazon:
        print(f"Ошибка: число {A} находится вне рабочего диапазона (0-{work_diapazon-1})")
    else:
        print(f"Число {A} находится в рабочем диапазоне (0-{work_diapazon-1})")
    
    return A

def check_with_opss(osn_inform, osn_kontrol, chisl):
    print("\nПроверка переводом числа в ОПСС")
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    osn_kolvo = len(osn)
    
    print(f"Основания системы: {osn}")
    print(f"Остатки числа: {chisl}")

    #1. Вычисление только необходимых констант c_ij
    print("Вычисление констант c_ij (обратные элементы)")
    c = numpy.zeros((osn_kolvo, osn_kolvo), dtype=object)
    
    #Вычисляем только элементы ниже главной диагонали
    for i in range(osn_kolvo):
        for j in range(i + 1, osn_kolvo):
            try:
                c[i][j] = pow(int(osn[i]), -1, int(osn[j]))
                print(f"c[{i}][{j}] = inv({osn[i]}) mod {osn[j]} = {c[i][j]}")
            except ValueError:
                print(f"Ошибка: невозможно найти c[{i}][{j}]")
                continue

    #2. Перевод в ОПСС с использованием вычисленных констант
    print("Перевод числа в ОПСС")
    opss = numpy.zeros(osn_kolvo, dtype=object)
    
    opss[0] = chisl[0] % osn[0]
    print(f"opss[0] = {chisl[0]} mod {osn[0]} = {opss[0]}")
    
    # Последующие разряды
    for i in range(1, osn_kolvo):
        temp = chisl[i]
        for j in range(i):
            temp = (temp - opss[j]) * c[j][i] % osn[i]
        opss[i] = temp % osn[i]
        print(f"opss[{i}] = {temp} mod {osn[i]} = {opss[i]}")
    
    print("\nЧисло в ОПСС:", opss)
    
    #Проверка контрольных остатков
    print("\nПроверка контрольных остатков:")
    error = False
    for i in range(len(osn_inform), osn_kolvo):
        if opss[i] != 0:
            print(f"Ошибка: opss[{i}] = {opss[i]} ≠ 0 (по контрольному основанию {osn[i]})")
            error = True
    if not error:
        print("Все контрольные остатки равны 0 - ошибок не обнаружено")
    
    #3. Перевод в десятичную систему
    print("Перевод из ОПСС в десятичную систему")
    opss_check = 0
    product = 1
    for i in range(osn_kolvo):
        term = opss[i] * product
        print(f"{opss[i]} * {product} = {term}")
        opss_check += term
        if i < osn_kolvo - 1:
            product *= osn[i]
    
    print(f"\nИтоговое число: {opss_check}")
    return opss_check

def check_with_combined_method(osn_inform, osn_kontrol, chisl):
    print("\nПроверка комбинированным методом ортогональных базисов и ОПСС")
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    osn_kolvo = len(osn)
    full_diapazon = numpy.prod(osn.astype(object))
    
    #1. Вычисление ортогональных базисов
    print("Вычисление ортогональных базисов")
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
        try:
            m[i] = pow(int(beta[i]), -1, int(osn[i]))
            print(f"m[{i}] = inv({beta[i]}) mod {osn[i]} = {m[i]}")
        except ValueError:
            print(f"Ошибка: невозможно найти m[{i}]")
            return None

    B = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        B[i] = m[i] * P[i]
        print(f"B[{i}] = {m[i]} * {P[i]} = {B[i]}")

    #2. Перевод базисов в СОК
    print("Перевод базисов в СОК")
    B_modular = numpy.zeros((osn_kolvo, osn_kolvo), dtype=object)
    for i in range(osn_kolvo):
        B_modular[i][:] = numpy.mod(B, osn[i])
        print(f"B_modular[{i}] = B mod {osn[i]} = {B_modular[i][:]}")

    #3. Вычисление констант для ОПСС
    print("Вычисление констант c_ij для ОПСС")
    c = numpy.zeros((osn_kolvo, osn_kolvo), dtype=object)
    for i in range(osn_kolvo):
        for j in range(i + 1, osn_kolvo):
            try:
                c[i][j] = pow(int(osn[i]), -1, int(osn[j]))
                print(f"c[{i}][{j}] = inv({osn[i]}) mod {osn[j]} = {c[i][j]}")
            except ValueError:
                print(f"Ошибка: невозможно найти c[{i}][{j}]")
                continue

    #4. Разложение базисов в ОПСС
    print("Разложение базисов в ОПСС")
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
    print("Расчет числа в ОПСС")
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
    print("Учет переносов")
    opss_total = numpy.zeros(osn_kolvo, dtype=object)
    perenos = 0
    for i in range(osn_kolvo):
        total = opss[i] + perenos
        opss_total[i] = total % osn[i]
        new_perenos = total // osn[i]
        print(f"opss_total[{i}] = ({opss[i]} + {perenos}) mod {osn[i]} = {opss_total[i]}")
        print(f"perenos = ({opss[i]} + {perenos}) // {osn[i]} = {new_perenos}")
        perenos = new_perenos
    
    print("\nЧисло в ОПСС с учетом переносов:", opss_total)
    
    #7. Перевод в десятичную систему
    print("Перевод в десятичную систему")
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

#Вывод введенных данных
print("Информационные основания:", osn_inform)
print("Контрольные основания:", osn_kontrol)
print("Остатки числа:", chisl)
print("Рабочий диапазон:", numpy.prod(osn_inform.astype(object)))
print("Полный диапазон:", numpy.prod(osn.astype(object)))

#Проверка всеми методами
A_orth = check_with_orthogonal_bases(osn_inform, osn_kontrol, chisl)
A_opss = check_with_opss(osn_inform, osn_kontrol, chisl)
A_combined = check_with_combined_method(osn_inform, osn_kontrol, chisl)

#Сравнение результатов
print(f"Метод ортогональных базисов: {A_orth}")
print(f"Метод ОПСС: {A_opss}")
print(f"Комбинированный метод: {A_combined}")

if A_orth == A_opss == A_combined:
    print("Все методы дали одинаковый результат")
else:
    print("Обнаружены расхождения в результатах!")

#Замер ресурсов
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
execution_time = time.time() - start_time
print(f"Пиковое использование памяти: {peak / 1024:.2f} KB")
print(f"Время выполнения: {execution_time:.6f} секунд")
