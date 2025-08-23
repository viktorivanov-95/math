import numpy
import time
import tracemalloc

#Оценка размера файла
def estimate_file_size(full_diapazon, work_diapazon, osn):
    avg_line_length = len(f"i = 99999, a = {numpy.array([999]*len(osn))}, A_orth = 99999, A_opss = 99999, status = 'OK'\n")
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

#Проверка чисел методом ортогональных базисов
def optimized_orthogonal_bases(osn_inform, osn_kontrol, chisl):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    osn_kolvo = len(osn)
    full_diapazon = numpy.prod(osn.astype(object))
    
    #1. Вычисление Pi = P / pi
    P = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        P[i] = full_diapazon // osn[i]
    
    #2. Вычисление βi = Pi mod pi
    beta = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        beta[i] = P[i] % osn[i]
    
    #3. Нахождение mi = βi^(-1) mod pi
    m = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        try:
            m[i] = pow(int(beta[i]), -1, int(osn[i]))
        except ValueError:
            return None
    
    #4. Вычисление базисов Bi = mi * Pi
    B = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        B[i] = m[i] * P[i]
    
    #5. Вычисление числа A
    A = 0
    for i in range(osn_kolvo):
        A = (A + chisl[i] * B[i]) % full_diapazon
    
    return A

#Проверка чисел методом перевода в ОПСС
def optimized_opss_conversion(osn_inform, osn_kontrol, chisl):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    osn_kolvo = len(osn)
    
    #1. Поиск констант c_ij (нижняя треугольная матрица)
    c = numpy.zeros((osn_kolvo-1, osn_kolvo), dtype=object)
    for j in range(osn_kolvo-1):
        for i in range(j+1):
            try:
                c[j][i] = pow(int(osn[i]), -1, int(osn[j+1]))
            except ValueError:
                return None
    
    #2. Расчет числа в ОПСС
    opss = numpy.zeros(osn_kolvo, dtype=object)
    opss[0] = chisl[0] % osn[0]
    
    if osn_kolvo > 1:
        opss[1] = ((chisl[1] - opss[0]) * c[0][0]) % osn[1]
    
    for i in range(2, osn_kolvo):
        opss[i] = (chisl[i] - opss[0]) * c[i-1][0]
        for j in range(1, i):
            opss[i] = (opss[i] - opss[j]) * c[i-1][j]
        opss[i] %= osn[i]
    
    #3. Перевод ОПСС в десятичное число
    opss_check = 0
    for i in range(osn_kolvo):
        opss_check += opss[i] * numpy.prod(osn[:i].astype(object)) if i > 0 else opss[i]
    
    return opss_check

#Генерация модулярных чисел
def generate_modular_numbers(osn_inform, osn_kontrol):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    full_diapazon = numpy.prod(osn.astype(object))
    numbers = []
    for i in range(int(full_diapazon)):
        a = numpy.mod(i, osn)
        numbers.append((i, a))
    return numbers

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
filename = 'modular_numbers_verification.txt'
save_to_file = input(f'Записать все числа с проверками в файл {filename}? [yes/no]: ').strip().lower()

file_created = False
results = []

if save_to_file == 'yes':
    file_created = True
    with open(filename, 'w', encoding='utf-8') as file:
        file.write(f'Информационные основания: {osn_inform}\n')
        file.write(f'Контрольные основания: {osn_kontrol}\n')
        file.write(f'Все основания: {osn}\n')
        file.write(f'Рабочий диапазон: 0-{work_diapazon-1}\n')
        file.write(f'Полный диапазон: 0-{full_diapazon-1}\n\n')
        
        numbers = generate_modular_numbers(osn_inform, osn_kontrol)
        
        file.write('Результаты проверки чисел\n')
        for i, a in numbers:
            A_orth = optimized_orthogonal_bases(osn_inform, osn_kontrol, a)
            A_opss = optimized_opss_conversion(osn_inform, osn_kontrol, a)
            
            status = 'OK' if i == A_orth == A_opss else 'ERROR'
            if A_orth is not None and A_orth > work_diapazon:
                status = 'OUT_OF_RANGE'
            
            result_line = f"i = {i}, a = {a}, A_orth = {A_orth}, A_opss = {A_opss}, status = '{status}'\n"
            file.write(result_line)
            results.append((i, a, A_orth, A_opss, status))
    
    print(f'Файл успешно создан: {filename}')

#Запрос на вывод на экран
if not file_created:
    numbers = generate_modular_numbers(osn_inform, osn_kontrol)
    for i, a in numbers:
        A_orth = optimized_orthogonal_bases(osn_inform, osn_kontrol, a)
        A_opss = optimized_opss_conversion(osn_inform, osn_kontrol, a)
        
        status = 'OK' if i == A_orth == A_opss else 'ERROR'
        if A_orth is not None and A_orth > work_diapazon:
            status = 'OUT_OF_RANGE'
        
        results.append((i, a, A_orth, A_opss, status))

user_input = input(f'\nВывести все числа с проверками на экран? (Их количество: {full_diapazon}) [yes/no]: ').strip().lower()

if user_input == 'yes':
    print('\nРезультаты проверки чисел')
    for result in results:
        i, a, A_orth, A_opss, status = result
        print(f"i = {i}, a = {a}, A_orth = {A_orth}, A_opss = {A_opss}, status = '{status}'")

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
