import numpy
import time
import tracemalloc
import sys
from itertools import combinations, product

#Проверка, принадлежит ли число разрещенному диапазону методом ортогональных базисов
def check_number_range_orthogonal_bases(osn_inform, osn_kontrol, chisl):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    osn_kolvo = len(osn)
    #Вычисляем диапазоны
    work_diapazon = numpy.prod(osn_inform.astype(object))
    full_diapazon = numpy.prod(osn.astype(object))
    #Вычисляем ортогональные базисы
    P = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        P[i] = full_diapazon // osn[i]
    beta = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        beta[i] = P[i] % osn[i]
    m = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        m[i] = pow(int(beta[i]), -1, int(osn[i]))
    B = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        B[i] = m[i] * P[i]
    #Находим число A в десятичной системе
    A = 0
    for i in range(osn_kolvo):
        A += chisl[i] * B[i]
    
    A = A % full_diapazon
    #Проверяем диапазон
    return A < work_diapazon, A, work_diapazon, full_diapazon

def analyze_single_number_error_correction(input_number, inform_codes, osn, osn_inform):
    n = len(osn)
    results = {}
    
    #Проверяем, находится ли число в информационном диапазоне
    is_in_inform_range, number_index, work_diapazon, full_diapazon = check_number_range_orthogonal_bases(
        osn_inform, numpy.array([], dtype=object), input_number
    )
    
    if not is_in_inform_range:
        print(f"Число {input_number} находится в ЗАПРЕЩЕННОМ диапазоне!")
        print(f"Десятичное значение: {number_index}")
        print(f"Рабочий диапазон: 0-{work_diapazon-1}")
        return None
    
    print(f"Число {input_number} находится в РАЗРЕШЕННОМ диапазоне (десятичное значение: {number_index})")
    
    #Анализируем ошибки для каждой кратности
    for error_count in range(1, n + 1):
        print(f"\n{'='*50}")
        print(f"АНАЛИЗ ОШИБОК КРАТНОСТИ {error_count}")
        print(f"{'='*50}")
        
        correctable_errors = []
        uncorrectable_errors = []
        total_errors = 0
        
        #Генерируем все комбинации позиций для error_count ошибок
        for error_positions in combinations(range(n), error_count):
            #Генерируем все возможные значения ошибок в этих позициях
            error_values = []
            for pos in error_positions:
                #Все возможные значения ошибки для данного модуля (кроме правильного)
                error_values.append([val for val in range(osn[pos]) if val != input_number[pos]])
            
            #Декартово произведение всех возможных значений ошибок
            for error_combination in product(*error_values):
                total_errors += 1
                
                #Создаем ошибочный код
                erroneous_code = input_number.copy()
                for i, pos in enumerate(error_positions):
                    erroneous_code[pos] = error_combination[i]
                
                #Проверяем, попадает ли ошибочный код в информационный диапазон
                is_in_inform_range = any(numpy.array_equal(erroneous_code, c) for c in inform_codes)
                
                error_info = {
                    'positions': error_positions,
                    'values': error_combination,
                    'erroneous_code': erroneous_code.copy(),
                    'correctable': not is_in_inform_range
                }
                
                if not is_in_inform_range:
                    correctable_errors.append(error_info)
                else:
                    uncorrectable_errors.append(error_info)
        
        #Выводим результаты
        print(f"\nНЕИСПРАВЛЯЕМЫЕ ошибки ({len(uncorrectable_errors)}):")
        for error in uncorrectable_errors:
            print(f"  Позиции: {error['positions']}, Ошибочные значения: {error['values']}")
            print(f"  Ошибочный код: {error['erroneous_code']}")
        
        print(f"\nИСПРАВЛЯЕМЫЕ ошибки ({len(correctable_errors)}):")
        for error in correctable_errors:
            print(f"  Позиции: {error['positions']}, Ошибочные значения: {error['values']}")
            print(f"  Ошибочный код: {error['erroneous_code']}")
        
        #Статистика
        if total_errors > 0:
            correctable_percent = (len(correctable_errors) / total_errors) * 100
            uncorrectable_percent = (len(uncorrectable_errors) / total_errors) * 100
            
            print(f"\nСТАТИСТИКА для ошибок кратности {error_count}:")
            print(f"Всего возможных ошибок: {total_errors}")
            print(f"Исправляемые ошибки: {len(correctable_errors)} ({correctable_percent:.2f}%)")
            print(f"Неисправляемые ошибки: {len(uncorrectable_errors)} ({uncorrectable_percent:.2f}%)")
        
        results[error_count] = {
            'total': total_errors,
            'correctable': len(correctable_errors),
            'uncorrectable': len(uncorrectable_errors),
            'correctable_percent': correctable_percent if total_errors > 0 else 0,
            'uncorrectable_percent': uncorrectable_percent if total_errors > 0 else 0
        }
    
    return results

#Генерируем все коды в информационном диапазоне
def generate_inform_codes(osn, osn_inform):
    work_diapazon = 1
    for base in osn_inform:
        work_diapazon *= base
    
    inform_codes = []
    for i in range(work_diapazon):
        code = numpy.array([i % base for base in osn], dtype=numpy.int64)
        inform_codes.append(code)
    
    return inform_codes

#Оценка приблизительного размера файла
def estimate_file_size(osn, osn_inform, input_number):
    n = len(osn)
        avg_line_length = 100
    
    total_size = 0
    header_size = 500
    
    for error_count in range(1, n + 1):
        error_combinations = 1
        for pos in range(error_count):
            error_combinations *= (n - pos) // (pos + 1)
        
        error_values = 1
        for i in range(error_count):
            avg_error_values = sum(osn) / len(osn) - 1
            error_values *= avg_error_values
        
        total_errors = error_combinations * error_values
        total_size += total_errors * avg_line_length
    
    return total_size + header_size

def convert_bytes(size):
    for x in ['bytes', 'KB', 'MB', 'GB']:
        if size < 1024.0:
            return "%3.1f %s" % (size, x)
        size /= 1024.0

def main():
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
    
    #Объединение оснований
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    
    #Ввод анализируемого числа
    print(f"\nВведите число для анализа (формат: через пробел, {len(osn)} чисел):")
    input_str = input().strip()
    input_number = numpy.array([int(x) for x in input_str.split()], dtype=numpy.int64)
    
    #Проверка корректности ввода
    if len(input_number) != len(osn):
        print(f"Ошибка: нужно ввести ровно {len(osn)} чисел!")
        return
    
    for i, val in enumerate(input_number):
        if val >= osn[i]:
            print(f"Ошибка: число в позиции {i} должно быть меньше {osn[i]}!")
            return
    
    #Проверка диапазона методом ортогональных базисов
    print("\nПроверка диапазона числа методом ортогональных базисов...")
    is_in_inform_range, number_index, work_diapazon, full_diapazon = check_number_range_orthogonal_bases(
        osn_inform, osn_kontrol, input_number
    )
    
    if not is_in_inform_range:
        print(f"ОШИБКА: Число {input_number} находится в ЗАПРЕЩЕННОМ диапазоне!")
        print(f"Десятичное значение: {number_index}")
        print(f"Рабочий диапазон: 0-{work_diapazon-1}")
        print("Анализ ошибок возможен только для чисел из разрешенного диапазона.")
        return
    
    print(f"Число {input_number} находится в РАЗРЕШЕННОМ диапазоне (десятичное значение: {number_index})")
    
    #Оценка размера файла
    estimated_size = estimate_file_size(osn, osn_inform, input_number)
    print(f"\nПриблизительный размер файла с результатами: {convert_bytes(estimated_size)}")
    
    #Запрос на запись в файл
    filename = f'error_analysis_{"_".join(map(str, input_number))}.txt'
    save_to_file = input(f'\nЗаписать анализ ошибок в файл {filename}? [yes/no]: ').strip().lower()
    
    if save_to_file != 'yes':
        print("Запись в файл отменена.")
        return
    
    #Генерируем информационные коды
    print("\nГенерация информационных кодов...")
    inform_codes = generate_inform_codes(osn, osn_inform)
    print(f"Сгенерировано {len(inform_codes)} информационных кодов")
    
    #Анализируем ошибки
    with open(filename, 'w', encoding='utf-8') as file:
        #Перенаправляем вывод в файл
        original_stdout = sys.stdout
        sys.stdout = file
        
        try:
            print("ДЕТАЛЬНЫЙ АНАЛИЗ ОШИБОК ДЛЯ ЧИСЛА")
            print("=" * 60)
            print(f"Введенное число: {input_number}")
            print(f"Десятичное значение: {number_index}")
            print(f"Основания: {osn}")
            print(f"Информационные основания: {osn_inform}")
            print(f"Контрольные основания: {osn_kontrol}")
            print(f"Рабочий диапазон: 0-{work_diapazon-1}")
            print(f"Полный диапазон: 0-{full_diapazon-1}")
            print("=" * 60)
            
            results = analyze_single_number_error_correction(input_number, inform_codes, osn, osn_inform)
            
            if results:
                print("\n" + "="*70)
                print("ОБЩАЯ СТАТИСТИКА ИСПРАВЛЕНИЯ ОШИБОК")
                print("="*70)
                
                for error_count in range(1, len(osn) + 1):
                    stats = results[error_count]
                    print(f"\nОШИБКИ КРАТНОСТИ {error_count}:")
                    print(f"  Всего возможных ошибок: {stats['total']}")
                    print(f"  Исправляемые ошибки: {stats['correctable']} ({stats['correctable_percent']:.2f}%)")
                    print(f"  Неисправляемые ошибки: {stats['uncorrectable']} ({stats['uncorrectable_percent']:.2f}%)")
                    print(f"  Эффективность обнаружения: {stats['correctable_percent']:.2f}%")
        
        finally:
            sys.stdout = original_stdout
    
    print(f"\nАнализ завершен! Результаты сохранены в файле: {filename}")
    
    #Выводим краткую статистику на экран
    if results:
        print("\nКРАТКАЯ СТАТИСТИКА:")
        for error_count in range(1, len(osn) + 1):
            stats = results[error_count]
            print(f"Ошибки кратности {error_count}: {stats['correctable_percent']:.1f}% обнаруживаемых")

if __name__ == "__main__":
    main()
