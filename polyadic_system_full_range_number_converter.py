import numpy
import time
import tracemalloc

# Оценка размера файла
def estimate_file_size(full_diapazon, osn):
    avg_line_length = len(f"i = 99999, a = {numpy.array([999]*len(osn))}, A_orth = 99999, A_opss = 99999, status = 'OK'\n")
    header_lines = 8
    avg_header_line = 50
    
    total_size = (full_diapazon * avg_line_length + 
                 header_lines * avg_header_line)
    
    return total_size

def convert_bytes(size):
    for x in ['bytes', 'KB', 'MB', 'GB']:
        if size < 1024.0:
            return "%3.1f %s" % (size, x)
        size /= 1024.0

# Проверка чисел методом ортогональных базисов
def optimized_orthogonal_bases(bases, residues):
    bases_count = len(bases)
    full_diapazon = numpy.prod(bases.astype(object))
    
    # 1. Вычисление Pi = P / pi
    P = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        P[i] = full_diapazon // bases[i]
    
    # 2. Вычисление βi = Pi mod pi
    beta = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        beta[i] = P[i] % bases[i]
    
    # 3. Нахождение mi = βi^(-1) mod pi
    m = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        try:
            m[i] = pow(int(beta[i]), -1, int(bases[i]))
        except ValueError:
            return None
    
    # 4. Вычисление базисов Bi = mi * Pi
    B = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        B[i] = m[i] * P[i]
    
    # 5. Вычисление числа A
    A = 0
    for i in range(bases_count):
        A = (A + residues[i] * B[i]) % full_diapazon
    
    return A

# Проверка чисел методом перевода в ОПСС
def optimized_opss_conversion(bases, residues):
    bases_count = len(bases)
    
    # 1. Поиск констант c_ij (нижняя треугольная матрица)
    c = numpy.zeros((bases_count-1, bases_count), dtype=object)
    for j in range(bases_count-1):
        for i in range(j+1):
            try:
                c[j][i] = pow(int(bases[i]), -1, int(bases[j+1]))
            except ValueError:
                return None
    
    # 2. Расчет числа в ОПСС
    opss = numpy.zeros(bases_count, dtype=object)
    opss[0] = residues[0] % bases[0]
    
    if bases_count > 1:
        opss[1] = ((residues[1] - opss[0]) * c[0][0]) % bases[1]
    
    for i in range(2, bases_count):
        opss[i] = (residues[i] - opss[0]) * c[i-1][0]
        for j in range(1, i):
            opss[i] = (opss[i] - opss[j]) * c[i-1][j]
        opss[i] %= bases[i]
    
    # 3. Перевод ОПСС в десятичное число
    opss_check = 0
    for i in range(bases_count):
        if i > 0:
            opss_check += opss[i] * numpy.prod(bases[:i].astype(object))
        else:
            opss_check += opss[i]
    
    return opss_check

# Генерация модулярных чисел
def generate_modular_numbers(bases):
    full_diapazon = numpy.prod(bases.astype(object))
    numbers = []
    for i in range(int(full_diapazon)):
        a = numpy.mod(i, bases)
        numbers.append((i, a))
    return numbers

# Начало программы
start_time = time.time()
tracemalloc.start()

print("Введите основания системы (все основания должны быть попарно взаимно простыми):")
bases_count = int(input('Введите количество оснований: '))
bases = numpy.zeros(bases_count, dtype=object)

for i in range(bases_count):
    bases[i] = int(input(f'Введите {i+1} основание: '))

full_diapazon = numpy.prod(bases.astype(object))

print(f'\nПолный диапазон: 0-{full_diapazon-1} ({full_diapazon} чисел)')

# Оценка размера файла
estimated_size = estimate_file_size(full_diapazon, bases)
print(f"\nПриблизительный размер файла с результатами: {convert_bytes(estimated_size)}")

# Запрос на запись в файл
filename = 'modular_numbers_verification.txt'
save_to_file = input(f'Записать все числа с проверками в файл {filename}? [yes/no]: ').strip().lower()

file_created = False
results = []

if save_to_file == 'yes':
    file_created = True
    with open(filename, 'w', encoding='utf-8') as file:
        file.write(f'Основания системы: {bases}\n')
        file.write(f'Полный диапазон: 0-{full_diapazon-1}\n\n')
        
        numbers = generate_modular_numbers(bases)
        
        file.write('Результаты проверки чисел\n')
        for i, a in numbers:
            A_orth = optimized_orthogonal_bases(bases, a)
            A_opss = optimized_opss_conversion(bases, a)
            
            status = 'OK' if i == A_orth == A_opss else 'ERROR'
            
            result_line = f"i = {i}, a = {a}, A_orth = {A_orth}, A_opss = {A_opss}, status = '{status}'\n"
            file.write(result_line)
            results.append((i, a, A_orth, A_opss, status))
    
    print(f'Файл успешно создан: {filename}')

# Запрос на вывод на экран
if not file_created:
    numbers = generate_modular_numbers(bases)
    for i, a in numbers:
        A_orth = optimized_orthogonal_bases(bases, a)
        A_opss = optimized_opss_conversion(bases, a)
        
        status = 'OK' if i == A_orth == A_opss else 'ERROR'
        
        results.append((i, a, A_orth, A_opss, status))

user_input = input(f'\nВывести все числа с проверками на экран? (Их количество: {full_diapazon}) [yes/no]: ').strip().lower()

if user_input == 'yes':
    print('\nРезультаты проверки чисел')
    for result in results:
        i, a, A_orth, A_opss, status = result
        print(f"i = {i}, a = {a}, A_orth = {A_orth}, A_opss = {A_opss}, status = '{status}'")

# Замер памяти и времени
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
