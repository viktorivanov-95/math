import numpy
import time
import tracemalloc
import matplotlib.pyplot as plt
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

# Функция для запуска проверки с определенной точностью и сбора разниц
def run_precision_check(bases, precision, numbers):
    differences = []
    
    for i, a in numbers:
        A_relative = relative_values_method(bases, a, precision)
        if A_relative is not None:
            diff = abs(float(A_relative) - i)
            differences.append(diff)
        else:
            differences.append(float('nan'))  # Если ошибка, добавляем NaN
    
    return differences

# Функция для вычисления статистики разниц
def calculate_difference_stats(differences):
    # Фильтруем NaN значения
    valid_differences = [d for d in differences if not numpy.isnan(d)]
    
    if not valid_differences:
        return float('nan'), float('nan'), float('nan')
    
    min_diff = min(valid_differences)
    max_diff = max(valid_differences)
    avg_diff = sum(valid_differences) / len(valid_differences)
    
    return min_diff, max_diff, avg_diff

# Основная функция программы
def main():
    start_time = time.time()
    tracemalloc.start()
    
    print("Введите основания системы (все основания должны быть попарно взаимно простыми):")
    
    # Ввод данных
    bases_count = int(input('Введите количество оснований: '))
    bases = numpy.zeros(bases_count, dtype=object)
    
    for i in range(bases_count):
        bases[i] = int(input(f'Введите {i+1} основание: '))
    
    # Ввод диапазона точности
    print("\nНастройка диапазона точности метода относительных величин")
    while True:
        try:
            min_precision = int(input('Введите минимальное количество знаков после запятой (1-50): '))
            max_precision = int(input('Введите максимальное количество знаков после запятой (1-50): '))
            if (1 <= min_precision <= 50 and 1 <= max_precision <= 50 and 
                min_precision <= max_precision):
                break
            else:
                print("Ошибка: введите числа от 1 до 50, минимальное ≤ максимальному")
        except ValueError:
            print("Ошибка: введите целые числа")
    
    # Генерация чисел в модулярной системе
    numbers, full_diapazon = generate_modular_numbers(bases)
    decimal_numbers = [i for i, a in numbers]  # Десятичные числа для оси X
    
    print(f'\nПолный диапазон: 0-{full_diapazon-1} ({full_diapazon} чисел)')
    
    # Сбор данных для разных значений точности
    precision_differences = {}
    precision_stats = {}
    
    print(f"\nЗапуск проверки для диапазона точности {min_precision}-{max_precision} знаков...")
    
    for precision in range(min_precision, max_precision + 1):
        print(f"Запуск проверки для точности {precision} знаков после запятой...")
        differences = run_precision_check(bases, precision, numbers)
        precision_differences[precision] = differences
        
        # Вычисляем статистику для текущей точности
        min_diff, max_diff, avg_diff = calculate_difference_stats(differences)
        precision_stats[precision] = {
            'min': min_diff,
            'max': max_diff,
            'avg': avg_diff
        }
        
        print(f"Завершено для точности {precision}")
    
    # Вывод статистики по точности
    print("Статистика разниц метода")
    
    for precision, stats in precision_stats.items():
        if not numpy.isnan(stats['min']):
            print(f"Для приближения {precision}: минимальная = {stats['min']:.10f}, максимальная = {stats['max']:.10f}, средняя = {stats['avg']:.10f}")
        else:
            print(f"Для приближения {precision}: нет валидных данных")
    
    # Построение графика
    print("\nПостроение графика...")
    plt.figure(figsize=(14, 8))
    
    for precision, differences in precision_differences.items():
        plt.plot(decimal_numbers, differences, label=f'Точность {precision} знаков', linewidth=1, alpha=0.7)
    
    plt.xlabel('Десятичное число')
    plt.ylabel('Разница метода относительных величин')
    plt.title('Зависимость точности метода относительных величин от количества знаков после запятой')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Сохранение графика
    plot_filename = 'relative_method_accuracy_plot.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"График сохранен как {plot_filename}")
    
    # Сохранение статистики в файл
    stats_filename = 'relative_method_stats.txt'
    with open(stats_filename, 'w', encoding='utf-8') as f:
        f.write(f"Основания системы: {bases}\n")
        f.write(f"Полный диапазон: 0-{full_diapazon-1}\n")
        f.write(f"Диапазон точности: {min_precision}-{max_precision} знаков\n\n")
        
        for precision, stats in precision_stats.items():
            if not numpy.isnan(stats['min']):
                f.write(f"Для приближения {precision}: минимальная = {stats['min']:.10f}, максимальная = {stats['max']:.10f}, средняя = {stats['avg']:.10f}\n")
            else:
                f.write(f"Для приближения {precision}: нет валидных данных\n")
    
    print(f"Статистика сохранена в файл: {stats_filename}")
    
    # Показать график
    plt.show()
    
    # Замер ресурсов
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    execution_time = time.time() - start_time
    
    print("МЕТРИКИ ВЫПОЛНЕНИЯ")
    print(f'Пиковое использование памяти: {peak / 1024:.2f} KB')
    print(f'Время выполнения: {execution_time:.6f} секунд')
    print(f'Диапазон точности: {min_precision}-{max_precision} знаков после запятой')

if __name__ == "__main__":
    main()
