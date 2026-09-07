import numpy
import time
import tracemalloc

# Вычисляем размер файла
def estimate_file_size(full_diapazon, osn):
    # Средняя длина строки
    avg_line = len(f"Проверка числа i = 99999, a = {numpy.array([999]*len(osn))} \
                Вычисленное десятичное число A = 99999")
    
    # Количество строк заголовка
    header_lines = 6
    avg_header_line = 50  # средняя длина строки заголовка
    
    total_size = (full_diapazon * avg_line + header_lines * avg_header_line)
    
    return total_size

def convert_bytes(size):
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024.0:
            return "%3.1f %s" % (size, x)
        size /= 1024.0

# Генерация модулярных чисел внутри заданного диапазона
def generate_modular_numbers(osn):
    full_diapazon = numpy.prod(osn.astype(object))
    numbers = []
    for i in range(int(full_diapazon)):
        a = numpy.mod(i, osn)
        numbers.append((i, a))
    return numbers, full_diapazon

def modular_to_decimal(bases, residues):
    bases_count = len(bases)
    full_diapazon = numpy.prod(bases.astype(object))
    
    # Вычисление ортогональных базисов
    P = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        P[i] = full_diapazon // bases[i]
    
    beta = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        beta[i] = P[i] % bases[i]
    
    m = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        try:
            m[i] = pow(int(beta[i]), -1, int(bases[i]))
        except ValueError:
            print(f"Ошибка: невозможно найти обратный элемент для beta[{i}] = {beta[i]} по модулю {bases[i]}")
            return None
    
    B = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        B[i] = m[i] * P[i]
    
    # Нахождение числа A
    A = 0
    for i in range(bases_count):
        A = A + residues[i] * B[i]
    A = A % full_diapazon
    
    return A

# Ввод данных
print("Введите основания системы (все основания должны быть попарно взаимно простыми):")
osn_kolvo = int(input('Введите количество оснований: '))
osn = numpy.zeros(osn_kolvo, dtype=object)

for i in range(osn_kolvo):
    osn[i] = int(input(f'Введите {i+1} основание: '))

# Генерация чисел в модулярной системе
numbers, full_diapazon = generate_modular_numbers(osn)

# Вывод информации о диапазоне
print(f'\nКоличество всех чисел (0-{full_diapazon-1}): {full_diapazon}')

# Оценка размера файла
estimated_size = estimate_file_size(full_diapazon, osn)
print(f"\nПриблизительный размер файла с результатами: {convert_bytes(estimated_size)}")

# Запрос на запись в файл
filename = 'modular_numbers_results.txt'
save_to_file = input(f'Записать все числа с результатами в файл {filename}? [yes/no]: ').strip().lower()

file_created = False

if save_to_file == 'yes':
    # Запись данных в файл
    with open(filename, 'w', encoding='utf-8') as file:
        file_created = True
        
        # Запись заголовочной информации
        file.write(f'Основания системы: {osn}\n')
        file.write(f'Полный диапазон: 0-{full_diapazon-1}\n\n')
        file.write(f'Количество всех чисел: {full_diapazon}\n\n')
        
        # Запись результатов перевода
        for i, a in numbers:
            A = modular_to_decimal(osn, a)
            if A is not None:
                file.write(f"i = {i}, a = {a}, A = {A}\n")
                if i == A:
                    file.write("Число корректно преобразовано из модулярной системы в десятичную.\n")
                else:
                    file.write("Ошибка преобразования: i и A не совпадают!\n")
    
    print(f'Файл успешно создан: {filename}')

# Запрос на вывод на экран
user_input = input(f'\nВывести все числа с результатами на экран? (Их количество: {full_diapazon}) [yes/no]: ').strip().lower()

if user_input == 'yes':
    for i, a in numbers:
        print(f"\nПроверка числа i = {i}, a = {a}:")
        A = modular_to_decimal(osn, a)
        if A is not None:
            print(f"Вычисленное десятичное число A = {A}")
            
            if i == A:
                print("Число корректно преобразовано из модулярной системы в десятичную.")
            else:
                print("Ошибка преобразования: i и A не совпадают!")

if file_created:
    print(f'\nПрограмма завершена. Результаты сохранены в файле: {filename}')
else:
    print('\nПрограмма завершена. Файл не создавался.')
