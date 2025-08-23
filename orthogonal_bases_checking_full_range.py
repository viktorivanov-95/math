import numpy
import time
import tracemalloc

#Вычисляем размер файла
def estimate_file_size(full_diapazon, work_diapazon, osn):
    #Средняя длина строки
    avg_line = len(f"Проверка числа i = 99999, a = {numpy.array([999]*len(osn))} \
                Вычисленное десятичное число A = 99999 \
                Результат проверки: a = {numpy.array([999]*len(osn))} - число не содержит ошибки (ошибок) \
                Число корректно преобразовано из модулярной системы в десятичную.")
    
    #Количество строк заголовка
    header_lines = 10
    avg_header_line = 50  #средняя длина строки заголовка
    
    total_size = (full_diapazon * avg_line + header_lines * avg_header_line)
    
    return total_size

def convert_bytes(size):
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024.0:
            return "%3.1f %s" % (size, x)
        size /= 1024.0

#Генерация модулярных чисел внутри заданного диапазона
def generate_modular_numbers(osn_inform, osn_kontrol):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    work_diapazon = numpy.prod(osn_inform.astype(object))
    full_diapazon = numpy.prod(osn.astype(object))
    numbers = []
    for i in range(int(full_diapazon)):
        a = numpy.mod(i, osn)
        numbers.append((i, a))
    return numbers, work_diapazon, full_diapazon

def check_number_with_orthogonal_bases(osn_inform, osn_kontrol, chisl):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    osn_kolvo = len(osn)
    work_diapazon = numpy.prod(osn_inform.astype(object))
    full_diapazon = numpy.prod(osn.astype(object))

    #Вычисление ортогональных базисов
    P = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        P[i] = full_diapazon // osn[i]

    beta = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        beta[i] = P[i] % osn[i]

    m = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        try:
            m[i] = pow(int(beta[i]), -1, int(osn[i]))
        except ValueError:
            print(f"Ошибка: невозможно найти обратный элемент для beta[{i}] = {beta[i]} по модулю {osn[i]}")
            return None, "Ошибка при вычислении обратного элемента"

    B = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        B[i] = m[i] * P[i]

    #Нахождение числа A
    A = 0
    for i in range(osn_kolvo):
        A = A + chisl[i] * B[i]
    A = A % full_diapazon

    #Проверка числа
    if A > work_diapazon:
        return A, f"{chisl} - представление числа, содержащее ошибку (ошибки)"
    else:
        return A, f"{chisl} - число не содержит ошибки (ошибок)"

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

#Генерация чисел в модулярной системе
numbers, work_diapazon, full_diapazon = generate_modular_numbers(osn_inform, osn_kontrol)
osn = numpy.concatenate((osn_inform, osn_kontrol))
#Вывод информации о диапазонах
print(f'\nКоличество чисел в разрешенном диапазоне (0-{work_diapazon-1}): {work_diapazon}')
print(f'Количество чисел в запрещенном диапазоне ({work_diapazon}-{full_diapazon-1}): {full_diapazon - work_diapazon}')
print(f'Количество всех чисел (0-{full_diapazon-1}): {full_diapazon}')
#Оценка размера файла
estimated_size = estimate_file_size(full_diapazon, work_diapazon, osn)
print(f"\nПриблизительный размер файла с результатами: {convert_bytes(estimated_size)}")
#Запрос на запись в файл
filename = 'modular_numbers_check_results.txt'
save_to_file = input(f'Записать все числа с результатами проверки в файл {filename}? [yes/no]: ').strip().lower()

file_created = False

if save_to_file == 'yes':
    #Запись данных в файл
    with open(filename, 'w', encoding='utf-8') as file:
        file_created = True
        
        #Запись заголовочной информации
        file.write(f'Информационные основания: {osn_inform}\n')
        file.write(f'Контрольные основания: {osn_kontrol}\n')
        file.write(f'Все основания: {osn}\n')
        file.write(f'Рабочий диапазон: 0-{work_diapazon-1}\n')
        file.write(f'Полный диапазон: 0-{full_diapazon-1}\n\n')
        file.write(f'Количество чисел в разрешенном диапазоне: {work_diapazon}\n')
        file.write(f'Количество чисел в запрещенном диапазоне: {full_diapazon - work_diapazon}\n')
        file.write(f'Количество всех чисел: {full_diapazon}\n\n')
        
        #Запись результатов проверки
        for i, a in numbers:
            A, result = check_number_with_orthogonal_bases(osn_inform, osn_kontrol, a)
            if A is not None:
                file.write(f"i = {i}, a = {a}, A = {A}, result = {result}\n")
                if i == A:
                    file.write("Число корректно преобразовано из модулярной системы в десятичную.\n")
                else:
                    file.write("Ошибка преобразования: i и A не совпадают!\n")
    
    print(f'Файл успешно создан: {filename}')

#Запрос на вывод на экран
user_input = input(f'\nВывести все числа с результатами проверки на экран? (Их количество: {full_diapazon}) [yes/no]: ').strip().lower()

if user_input == 'yes':
    print('\n===== Результаты проверки чисел =====')
    for i, a in numbers:
        print(f"\nПроверка числа i = {i}, a = {a}:")
        A, result = check_number_with_orthogonal_bases(osn_inform, osn_kontrol, a)
        if A is not None:
            print(f"Вычисленное десятичное число A = {A}")
            print(f"Результат проверки: {result}")
            
            if i == A:
                print("Число корректно преобразовано из модулярной системы в десятичную.")
            else:
                print("Ошибка преобразования: i и A не совпадают!")

if file_created:
    print(f'\nПрограмма завершена. Результаты сохранены в файле: {filename}')
else:
    print('\nПрограмма завершена. Файл не создавался.')
