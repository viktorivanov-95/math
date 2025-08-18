import numpy
import math

#Вычисляем минимальные/максимальные расстояния между кодовыми словами
def calculate_distances(code, codes_to_compare):
    min_dist = float('inf')
    max_dist = 0

    #Пропускаем, чтобы не сравнивать одинаковые кодовые слова
    for other_code in codes_to_compare:
        if numpy.array_equal(code, other_code):
            continue
        
        dist = numpy.sum(code != other_code)

        if dist < min_dist:
            min_dist = dist
        if dist > max_dist:
            max_dist = dist
    
    return min_dist, max_dist


#Вычисляем размер файла
def estimate_file_size(full_diapazon, work_diapazon, osn):
    def get_avg_i_length(start, end):
        #Если максимальное значение числа i <= 10^10, вычисляем
        #среднее кол-во цифр точно. В противном случае, через логарифм,
        #т.к. логарифм вычисляется быстрее
        if end <= 10**10:
            return (len(str(start)) + len(str(end))) / 2
        else:
            return math.log10(end) + 1
        
    min_a_length = len(f"a = {numpy.array([0] * len(osn))}")
    max_a_length = len(f"a = {numpy.array([x-1 for x in osn])}")
    avg_a_length = (min_a_length + max_a_length) / 2
    max_possible_dist = len(osn)
    avg_dist_length = len(f"min_dist = 0, max_dist = {max_possible_dist}")
    
    def avg_line_length(start, end):
        avg_i_len = get_avg_i_length(start, end)
        return (
            len("i = ") 
            + avg_i_len 
            + len(", ") 
            + avg_a_length 
            + len(", ") 
            + avg_dist_length 
            + 1  # для "\n"
        )

    avg_line_work = avg_line_length(0, work_diapazon - 1)
    avg_line_control = avg_line_length(work_diapazon, full_diapazon - 1)
    header_size = 15 * 50
    total_size = (
        work_diapazon * avg_line_work 
        + (full_diapazon - work_diapazon) * avg_line_control 
        + header_size
    )
    return total_size

def convert_bytes(size):
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024:
            return "%3.1f %s" % (size, x)
        size /= 1024

#Ввод данных с использованием Python int (бесконечная точность)
osn_inform_kolvo = int(input('Введите количество информационных оснований: '))
osn_kontrol_kolvo = int(input('Введите количество контрольных оснований: '))
osn_inform = []
osn_kontrol = []

for i in range(osn_inform_kolvo):
    osn_inform.append(int(input(f'Введите {i+1} информационное основание: ')))
for i in range(osn_kontrol_kolvo):
    osn_kontrol.append(int(input(f'Введите {i+1} контрольное основание: ')))

#Объединение оснований (используем Python списки для избежания переполнения)
osn = osn_inform + osn_kontrol

#Вычисление диапазонов с помощью math.prod (без переполнения)
work_diapazon = math.prod(osn_inform)
full_diapazon = math.prod(osn)

print(f'\nКоличество чисел в разрешенном диапазоне (0-{work_diapazon-1}): {work_diapazon}')
print(f'Количество чисел в запрещенном диапазоне ({work_diapazon}-{full_diapazon-1}): {full_diapazon - work_diapazon}')
print(f'Количество всех чисел (0-{full_diapazon-1}): {full_diapazon}')

osn_numpy = numpy.array(osn, dtype=object)  # Для numpy.mod

inform_min_dists = []
inform_max_dists = []
control_min_dists = []
control_max_dists = []

#Оценка размера файла
estimated_size = estimate_file_size(full_diapazon, work_diapazon, osn)
print(f"\nПриблизительный размер файла с результатами: {convert_bytes(estimated_size)}")

#Запрос на запись в файл
save_to_file = input('Записать все числа с расстояниями в файл? [yes/no]: ').strip().lower()

filename = 'modular_numbers_output.txt'
file_created = False

if save_to_file == 'yes':
    #Запись данных в файл
    with open(filename, 'w', encoding='utf-8') as file:
        file_created = True
        
        file.write(f'Информационные основания: {osn_inform}\n')
        file.write(f'Контрольные основания: {osn_kontrol}\n')
        file.write(f'Все основания: {osn}\n')
        file.write(f'Рабочий диапазон: 0-{work_diapazon-1}\n')
        file.write(f'Полный диапазон: 0-{full_diapazon-1}\n\n')
        file.write(f'Количество чисел в разрешенном диапазоне: {work_diapazon}\n')
        file.write(f'Количество чисел в запрещенном диапазоне: {full_diapazon - work_diapazon}\n')
        file.write(f'Количество всех чисел: {full_diapazon}\n\n')
        
        overall_inform_min = float('inf')
        overall_inform_max = 0
        overall_control_min = float('inf')
        overall_control_max = 0

        inform_codes = []
        for i in range(work_diapazon):
            inform_codes.append(numpy.mod(i, osn_numpy))
        
        #Запись разрешенного диапазона с расстояниями
        file.write('\n===== Разрешенный диапазон (информационный) с расстояниями =====\n')
        for i in range(work_diapazon):
            code = numpy.mod(i, osn)
            min_dist, max_dist = calculate_distances(code, inform_codes)
            inform_min_dists.append(min_dist)
            inform_max_dists.append(max_dist)
            file.write(f'i = {i}, a = {code}, min_dist = {min_dist}, max_dist = {max_dist}\n')
            
            if min_dist < overall_inform_min:
                overall_inform_min = min_dist
            if max_dist > overall_inform_max:
                overall_inform_max = max_dist
        
        #Запись запрещенного диапазона с расстояниями
        file.write('\n===== Запрещенный диапазон (контрольный) с расстояниями до информационных кодов =====\n')
        for i in range(work_diapazon, full_diapazon):
            code = numpy.mod(i, osn)
            min_dist, max_dist = calculate_distances(code, inform_codes)
            control_min_dists.append(min_dist)
            control_max_dists.append(max_dist)
            file.write(f'i = {i}, a = {code}, min_dist = {min_dist}, max_dist = {max_dist}\n')
            
            if min_dist < overall_control_min:
                overall_control_min = min_dist
            if max_dist > overall_control_max:
                overall_control_max = max_dist
        
        file.write('\n===== Общая статистика по расстояниям =====\n')
        file.write(f'Минимальное расстояние в разрешенном диапазоне: {overall_inform_min}\n')
        file.write(f'Максимальное расстояние в разрешенном диапазоне: {overall_inform_max}\n')
        file.write(f'Минимальное расстояние от запрещенного диапазона до информационных кодов: {overall_control_min}\n')
        file.write(f'Максимальное расстояние от запрещенного диапазона до информационных кодов: {overall_control_max}\n')
    
    print(f'Файл успешно создан: {filename}')
else:
    #Если файл не создается, все равно нужно вычислить расстояния для вывода на экран
    inform_codes = []
    for i in range(work_diapazon):
        inform_codes.append(numpy.mod(i, osn_numpy))

    for i in range(work_diapazon):
        code = numpy.mod(i, osn)
        min_dist, max_dist = calculate_distances(code, inform_codes)
        inform_min_dists.append(min_dist)
        inform_max_dists.append(max_dist)
    
    for i in range(work_diapazon, full_diapazon):
        code = numpy.mod(i, osn)
        min_dist, max_dist = calculate_distances(code, inform_codes)
        control_min_dists.append(min_dist)
        control_max_dists.append(max_dist)
    
    print('Запись в файл пропущена.')

#Запрос на вывод на экран
user_input = input(f'\nВывести все числа с расстояниями на экран? (Их количество: {full_diapazon}) [yes/no]: ').strip().lower()

if user_input == 'yes':
    print('\n===== Статистика по расстояниям =====')
    print(f'Минимальное расстояние в разрешенном диапазоне: {min(inform_min_dists)}')
    print(f'Максимальное расстояние в разрешенном диапазоне: {max(inform_max_dists)}')
    print(f'Минимальное расстояние от запрещенного диапазона до информационных кодов: {min(control_min_dists)}')
    print(f'Максимальное расстояние от запрещенного диапазона до информационных кодов: {max(control_max_dists)}')
    
    print('\n===== Разрешенный диапазон (информационный) с расстояниями =====')
    for i in range(work_diapazon):
        code = numpy.mod(i, osn)
        min_dist, max_dist = inform_min_dists[i], inform_max_dists[i]
        print(f'i = {i}, a = {code}, min_dist = {min_dist}, max_dist = {max_dist}')
    
    print('\n===== Запрещенный диапазон (контрольный) с расстояниями до информационных кодов =====')
    for i in range(work_diapazon, full_diapazon):
        code = numpy.mod(i, osn)
        idx = i - work_diapazon
        min_dist, max_dist = control_min_dists[idx], control_max_dists[idx]
        print(f'i = {i}, a = {code}, min_dist = {min_dist}, max_dist = {max_dist}')
else:
    print('Вывод чисел на экран пропущен.')

if file_created:
    print(f'\nПрограмма завершена. Результаты сохранены в файле: {filename}')
else:
    print('\nПрограмма завершена. Файл не создавался.')
