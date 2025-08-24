import numpy
import time
import tracemalloc

#Генерируем модулярные числа
def generate_modular_numbers(osn_inform, osn_kontrol):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    work_diapazon = numpy.prod(osn_inform.astype(object))
    full_diapazon = numpy.prod(osn.astype(object))
    numbers = []
    for i in range(int(full_diapazon)):
        a = numpy.mod(i, osn)
        numbers.append((i, a))
    return numbers, work_diapazon, full_diapazon

#Проверка методом ортогональных базисов
def check_with_orthogonal_bases(osn_inform, osn_kontrol, chisl):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    osn_kolvo = len(osn)
    work_diapazon = numpy.prod(osn_inform.astype(object))
    full_diapazon = numpy.prod(osn.astype(object))

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
            return None

    B = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        B[i] = m[i] * P[i]

    A = 0
    for i in range(osn_kolvo):
        A = A + chisl[i] * B[i]
    A = A % full_diapazon

    return A

#Проверка методом перевода в ОПСС
def check_with_opss(osn_inform, osn_kontrol, chisl):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    osn_kolvo = len(osn)
    
    #Вычисление констант c_ij (обратные элементы)
    c = numpy.zeros((osn_kolvo, osn_kolvo), dtype=object)
    for i in range(osn_kolvo):
        for j in range(i + 1, osn_kolvo):
            try:
                c[i][j] = pow(int(osn[i]), -1, int(osn[j]))
            except ValueError:
                continue
    
    #Перевод в ОПСС
    opss = numpy.zeros(osn_kolvo, dtype=object)
    opss[0] = chisl[0] % osn[0]
    
    #Последующие разряды
    for i in range(1, osn_kolvo):
        temp = chisl[i]
        for j in range(i):
            temp = (temp - opss[j]) * c[j][i] % osn[i]
        opss[i] = temp % osn[i]
    
    #Перевод из ОПСС в десятичную систему
    opss_check = 0
    product = 1
    for i in range(osn_kolvo):
        opss_check += opss[i] * product
        if i < osn_kolvo - 1:
            product *= osn[i]
    
    return opss_check

#Проверка методом совместного использования ортогональных базисов и ОПСС
def check_with_combined_method(osn_inform, osn_kontrol, chisl):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    osn_kolvo = len(osn)
    full_diapazon = numpy.prod(osn.astype(object))
    
    #1. Вычисление ортогональных базисов
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
            return None

    B = numpy.zeros(osn_kolvo, dtype=object)
    for i in range(osn_kolvo):
        B[i] = m[i] * P[i]
    
    #2. Перевод базисов в СОК
    B_modular = numpy.zeros((osn_kolvo, osn_kolvo), dtype=object)
    for i in range(osn_kolvo):
        B_modular[i][:] = B % osn[i]
    
    #3. Вычисление констант для ОПСС
    c = numpy.zeros((osn_kolvo, osn_kolvo), dtype=object)
    for i in range(osn_kolvo):
        for j in range(i + 1, osn_kolvo):
            try:
                c[i][j] = pow(int(osn[i]), -1, int(osn[j]))
            except ValueError:
                continue
    
    #4. Разложение базисов в ОПСС
    B_opss = numpy.zeros((osn_kolvo, osn_kolvo), dtype=object)
    for k in range(osn_kolvo):
        B_opss[k][0] = B_modular[k][0] % osn[0]
        
        for i in range(1, osn_kolvo):
            temp = B_modular[k][i]
            for j in range(i):
                temp = (temp - B_opss[k][j]) * c[j][i] % osn[i]
            B_opss[k][i] = temp % osn[i]

    #5. Расчет числа в ОПСС
    B_opss_chisl = numpy.zeros((osn_kolvo, osn_kolvo), dtype=object)
    for i in range(osn_kolvo):
        for j in range(osn_kolvo):
            B_opss_chisl[i][j] = chisl[i] * B_opss[i][j]
    
    opss = numpy.zeros(osn_kolvo, dtype=object)
    for j in range(osn_kolvo):
        for i in range(osn_kolvo):
            opss[j] += B_opss_chisl[i][j]
    
    #6. Учет переносов
    opss_total = numpy.zeros(osn_kolvo, dtype=object)
    perenos = 0
    for i in range(osn_kolvo):
        total = opss[i] + perenos
        opss_total[i] = total % osn[i]
        perenos = total // osn[i]
    
    #7. Перевод в десятичную систему
    A_combined = 0
    product = 1
    for i in range(osn_kolvo):
        A_combined += opss_total[i] * product
        if i < osn_kolvo - 1:
            product *= osn[i]
    
    return A_combined

#Оценка размера файла
def estimate_file_size(full_diapazon):
    """Оценивает приблизительный размер файла"""
    avg_line_length = 150  # Средняя длина строки с результатами проверки
    header_lines = 10      # Количество строк заголовка
    avg_header_length = 50 # Средняя длина строки заголовка
    
    total_size = (full_diapazon * avg_line_length + 
                 header_lines * avg_header_length)
    
    return total_size

def convert_bytes(size):
    """Конвертирует байты в удобочитаемый формат"""
    for x in ['bytes', 'KB', 'MB', 'GB']:
        if size < 1024.0:
            return "%3.1f %s" % (size, x)
        size /= 1024.0

#Ввод данных
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

#Оценка размера файла
estimated_size = estimate_file_size(full_diapazon)
print(f"\nПриблизительный размер файла с результатами: {convert_bytes(estimated_size)}")

#Запрос на запись в файл
filename = 'modular_numbers_verification.txt'
save_to_file = input(f'Записать все числа с результатами проверки в файл {filename}? [yes/no]: ').strip().lower()

file_created = False

if save_to_file == 'yes':
    with open(filename, 'w', encoding='utf-8') as file:
        file_created = True
        
        file.write(f'Информационные основания: {osn_inform}\n')
        file.write(f'Контрольные основания: {osn_kontrol}\n')
        file.write(f'Рабочий диапазон: 0-{work_diapazon-1}\n')
        file.write(f'Полный диапазон: 0-{full_diapazon-1}\n\n')
        
        #Проверка каждого числа и запись в файл
        for i, a in numbers:
            file.write(f"\nПроверка числа i = {i}, a = {a}:\n")
            
            #Проверка методом ортогональных базисов
            A_orth = check_with_orthogonal_bases(osn_inform, osn_kontrol, a)
            file.write(f"Метод ортогональных базисов: A = {A_orth}\n")
            
            #Проверка методом ОПСС
            A_opss = check_with_opss(osn_inform, osn_kontrol, a)
            file.write(f"Метод ОПСС: A = {A_opss}\n")
            
            #Проверка комбинированным методом
            A_combined = check_with_combined_method(osn_inform, osn_kontrol, a)
            file.write(f"Комбинированный метод: A = {A_combined}\n")
            
            #Сравнение результатов
            if i == A_orth == A_opss == A_combined:
                file.write("Все методы дали корректный результат\n")
            else:
                file.write("Обнаружено несоответствие:\n")
                if i != A_orth:
                    file.write(f"  - Метод ортогональных базисов: ожидалось {i}, получено {A_orth}\n")
                if i != A_opss:
                    file.write(f"  - Метод ОПСС: ожидалось {i}, получено {A_opss}\n")
                if i != A_combined:
                    file.write(f"  - Комбинированный метод: ожидалось {i}, получено {A_combined}\n")
            
            #Проверка на ошибки (по рабочему диапазону)
            if A_orth is not None and A_orth > work_diapazon:
                file.write("Обнаружена ошибка (число вне рабочего диапазона)\n")
    
    print(f'Файл успешно создан: {filename}')
else:
    print('Запись в файл пропущена.')

#Запрос на вывод на экран
user_input = input(f'\nВывести все числа с результатами проверки на экран? (Их количество: {full_diapazon}) [yes/no]: ').strip().lower()

if user_input == 'yes':
    #Проверка каждого числа и вывод на экран
    for i, a in numbers:
        print(f"\nПроверка числа i = {i}, a = {a}:")
        
        A_orth = check_with_orthogonal_bases(osn_inform, osn_kontrol, a)
        print(f"Метод ортогональных базисов: A = {A_orth}")
        
        A_opss = check_with_opss(osn_inform, osn_kontrol, a)
        print(f"Метод ОПСС: A = {A_opss}")
        
        A_combined = check_with_combined_method(osn_inform, osn_kontrol, a)
        print(f"Комбинированный метод: A = {A_combined}")
        
        if i == A_orth == A_opss == A_combined:
            print("Все методы дали корректный результат")
        else:
            print("Обнаружено несоответствие:")
            if i != A_orth:
                print(f"  - Метод ортогональных базисов: ожидалось {i}, получено {A_orth}")
            if i != A_opss:
                print(f"  - Метод ОПСС: ожидалось {i}, получено {A_opss}")
            if i != A_combined:
                print(f"  - Комбинированный метод: ожидалось {i}, получено {A_combined}")
        
        if A_orth is not None and A_orth > work_diapazon:
            print("Обнаружена ошибка (число вне рабочего диапазона)")
else:
    print('Вывод чисел на экран пропущен.')

if file_created:
    print(f'\nПрограмма завершена. Результаты сохранены в файле: {filename}')
else:
    print('\nПрограмма завершена. Файл не создавался.')
