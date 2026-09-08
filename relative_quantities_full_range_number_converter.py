import numpy
import time
import tracemalloc
from decimal import Decimal, getcontext, ROUND_HALF_UP

# Генерируем модулярные числа
def generate_modular_numbers(bases):
    full_diapazon = numpy.prod(bases.astype(object))
    numbers = []
    for i in range(int(full_diapazon)):
        a = numpy.mod(i, bases)
        numbers.append((i, a))
    return numbers, full_diapazon

# Проверка методом ортогональных базисов
def check_with_orthogonal_bases(bases, residues):
    bases_count = len(bases)
    full_diapazon = numpy.prod(bases.astype(object))
    
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
    
    A = 0
    for i in range(bases_count):
        A = A + residues[i] * B[i]
    A = A % full_diapazon
    
    return A

# Проверка методом перевода в ОПСС
def check_with_opss(bases, residues):
    bases_count = len(bases)
    
    # Вычисление констант c_ij (обратные элементы)
    c = numpy.zeros((bases_count, bases_count), dtype=object)
    for i in range(bases_count):
        for j in range(i + 1, bases_count):
            try:
                c[i][j] = pow(int(bases[i]), -1, int(bases[j]))
            except ValueError:
                continue
    
    # Перевод в ОПСС
    opss = numpy.zeros(bases_count, dtype=object)
    opss[0] = residues[0] % bases[0]
    
    # Последующие разряды
    for i in range(1, bases_count):
        temp = residues[i]
        for j in range(i):
            temp = (temp - opss[j]) * c[j][i] % bases[i]
        opss[i] = temp % bases[i]
    
    # Перевод из ОПСС в десятичную систему
    opss_check = 0
    product = 1
    for i in range(bases_count):
        opss_check += opss[i] * product
        if i < bases_count - 1:
            product *= bases[i]
    
    return opss_check

# Метод относительных величин с поддержкой больших чисел
def relative_values_method(bases, residues, precision=10):
    bases_count = len(bases)
    
    # Используем Decimal для вычисления диапазонов
    full_diapazon = Decimal(1)
    for base in bases:
        full_diapazon *= Decimal(str(base))
    
    # Настройка точности Decimal
    getcontext().prec = precision + 50
    getcontext().rounding = ROUND_HALF_UP
    
    # 1. Вычисление Pi = P/pi
    P = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        base_decimal = Decimal(str(bases[i]))
        P_decimal = full_diapazon / base_decimal
        # Ограничиваем количество знаков после запятой
        P_decimal = P_decimal.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
        P[i] = P_decimal
    
    # 2. Вычисление коэффициентов k[i]
    k = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        # Вычисляем обратный элемент с целыми числами
        Pi_int = int(Decimal(P[i]).to_integral_value(rounding=ROUND_HALF_UP))
        base_int = int(bases[i])
        inv_Pi = pow(Pi_int, -1, base_int)
        # Вычисляем k[i] с ограниченной точностью
        k_decimal = Decimal(str(inv_Pi)) / Decimal(str(base_int))
        k_decimal = k_decimal.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
        k[i] = k_decimal
    
    # 3. Вычисление A/P
    A_relative = Decimal('0')
    for i in range(bases_count):
        term = Decimal(str(residues[i])) * k[i]
        term = term.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
        A_relative += term
        A_relative = A_relative.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
    
    # 4. Вычисление A
    fractional_part = A_relative % Decimal('1')
    fractional_part = fractional_part.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
    A_relative = fractional_part * full_diapazon
    A_relative = A_relative.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
    
    return A_relative

# Оценка размера файла
def estimate_file_size(full_diapazon):
    """Оценивает приблизительный размер файла"""
    avg_line_length = 180  # Для 3 методов
    header_lines = 8
    avg_header_length = 50
    
    total_size = (full_diapazon * avg_line_length + 
                 header_lines * avg_header_length)
    
    return total_size

def convert_bytes(size):
    """Конвертирует байты в удобочитаемый формат"""
    for x in ['bytes', 'KB', 'MB', 'GB']:
        if size < 1024.0:
            return "%3.1f %s" % (size, x)
        size /= 1024.0

# Основная функция программы
start_time = time.time()
tracemalloc.start()

print("Введите основания системы (все основания должны быть попарно взаимно простыми):")

# Ввод данных
bases_count = int(input('Введите количество оснований: '))
bases = numpy.zeros(bases_count, dtype=object)

for i in range(bases_count):
    bases[i] = int(input(f'Введите {i+1} основание: '))

# Ввод точности для метода относительных величин
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

# Генерация чисел в модулярной системе
numbers, full_diapazon = generate_modular_numbers(bases)

print(f'\nПолный диапазон: 0-{full_diapazon-1} ({full_diapazon} чисел)')

# Оценка размера файла
estimated_size = estimate_file_size(full_diapazon)
print(f"\nПриблизительный размер файла с результатами: {convert_bytes(estimated_size)}")

# Запрос на запись в файл
filename = 'modular_numbers_verification.txt'
save_to_file = input(f'Записать все числа с результатами проверки в файл {filename}? [yes/no]: ').strip().lower()

file_created = False

if save_to_file == 'yes':
    with open(filename, 'w', encoding='utf-8') as file:
        file_created = True
        
        file.write(f'Основания системы: {bases}\n')
        file.write(f'Точность метода относительных величин: {precision} знаков после запятой\n')
        file.write(f'Полный диапазон: 0-{full_diapazon-1}\n\n')
        
        # Проверка каждого числа и запись в файл
        for i, a in numbers:
            file.write(f"\nПроверка числа i = {i}, a = {a}:\n")
            
            # Проверка методом ортогональных базисов
            A_orth = check_with_orthogonal_bases(bases, a)
            file.write(f"Метод ортогональных базисов: A = {A_orth}\n")
            
            # Проверка методом ОПСС
            A_opss = check_with_opss(bases, a)
            file.write(f"Метод ОПСС: A = {A_opss}\n")
            
            # Проверка методом относительных величин
            A_relative = relative_values_method(bases, a, precision)
            file.write(f"Метод относительных величин: A = {A_relative}\n")
            
            # Сравнение результатов точных методов
            if i == A_orth == A_opss:
                file.write("Точные методы дали корректный результат\n")
            else:
                file.write("Обнаружено несоответствие в точных методах:\n")
                if i != A_orth:
                    file.write(f"  - Метод ортогональных базисов: ожидалось {i}, получено {A_orth}\n")
                if i != A_opss:
                    file.write(f"  - Метод ОПСС: ожидалось {i}, получено {A_opss}\n")
            
            # Проверка точности метода относительных величин
            if A_relative is not None:
                diff = abs(float(A_relative) - i)
                file.write(f"Разница метода относительных величин: {diff:.10f}\n")
    
    print(f'Файл успешно создан: {filename}')
else:
    print('Запись в файл пропущена.')

# Запрос на вывод на экран
user_input = input(f'\nВывести все числа с результатами проверки на экран? (Их количество: {full_diapazon}) [yes/no]: ').strip().lower()

if user_input == 'yes':
    # Проверка каждого числа и вывод на экран
    for i, a in numbers:
        print(f"\nПроверка числа i = {i}, a = {a}:")
        
        A_orth = check_with_orthogonal_bases(bases, a)
        print(f"Метод ортогональных базисов: A = {A_orth}")
        
        A_opss = check_with_opss(bases, a)
        print(f"Метод ОПСС: A = {A_opss}")
        
        A_relative = relative_values_method(bases, a, precision)
        print(f"Метод относительных величин: A = {A_relative}")
        
        if i == A_orth == A_opss:
            print("Точные методы дали корректный результат")
        else:
            print("Обнаружено несоответствие в точных методах:")
            if i != A_orth:
                print(f"  - Метод ортогональных базисов: ожидалось {i}, получено {A_orth}")
            if i != A_opss:
                print(f"  - Метод ОПСС: ожидалось {i}, получено {A_opss}")
        
        if A_relative is not None:
            diff = abs(float(A_relative) - i)
            print(f"Разница метода относительных величин: {diff:.10f}")
else:
    print('Вывод чисел на экран пропущен.')

# Замер ресурсов
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
execution_time = time.time() - start_time

print(f'Пиковое использование памяти: {peak / 1024:.2f} KB')
print(f'Время выполнения: {execution_time:.6f} секунд')
print(f'Точность вычислений: {precision} знаков после запятой')

if file_created:
    print(f'Программа завершена. Результаты сохранены в файле: {filename}')
else:
    print('Программа завершена. Файл не создавался.')
