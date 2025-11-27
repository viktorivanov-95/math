import numpy
import math
from itertools import combinations, count, chain
from math import ceil


#Вычисление ортогональных базисов
def calculate_orthogonal_bases(osn):
    full_diapazon = math.prod([int(x) for x in osn])
    P = [full_diapazon // int(x) for x in osn]
    beta = [p % m for p, m in zip(P, osn)]
    m = [pow(int(b), -1, int(m)) if b != 0 else 0 for b, m in zip(beta, osn)]
    return [mi * pi for mi, pi in zip(m, P)]

#Проверка числа методом ортогональных базисов
def check_number(osn, number):
    B = calculate_orthogonal_bases(osn)
    total = sum(int(n) * int(b) for n, b in zip(number, B))
    A = total % math.prod([int(x) for x in osn])
    return A

#Вычисление обнаружения ошибки новым методом:
def new_method_decoding(osn_inform, osn_kontrol, chisl, max_error_correction):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    k = len(osn_kontrol)
    n = len(osn_inform)
    n_plus_k = len(osn)
    full_diapazon = math.prod([int(x) for x in osn])
    inform_diapazon = math.prod([int(x) for x in osn_inform])
    osn_inform_mult = math.prod([int(x) for x in osn_inform])
    osn_kontrol_mult = math.prod([int(x) for x in osn_kontrol])
    A = check_number(osn, chisl)
    if A < inform_diapazon:
        return chisl.tolist()
    else:
        #Проверка на ошибки только по информационным модулям
        inform_errors = A - int((A/osn_kontrol_mult))*osn_kontrol_mult
        inform_error_vector = [0 for x in range(n_plus_k)]
        for i in range(n_plus_k):
            inform_error_vector[i] = (int((A/osn_kontrol_mult))*osn_kontrol_mult)%osn[i]
        inform_error_vector_checking = sum(1 for x in inform_error_vector if x != 0)
        if inform_errors < osn_inform_mult and inform_error_vector_checking <= max_error_correction:
            print('Ошибка по информационным основаниям')
            #Вычитаем из самого числа вектор ошибки:
            print('Вектор ошибки = ', inform_error_vector)
            right_modular_vector = [(chisl[i] - inform_error_vector[i]) % osn[i] for i in range(len(osn))]
            print('Правильное число = ', right_modular_vector)
            print('Проверим его на правильность')
            A = check_number(osn, right_modular_vector)
            print('Число в позиционной системе равно = ', A)
            if A < osn_inform_mult:
                print('Число находится в разрешенном диапазоне, ошибка исправлена правильно')
                return right_modular_vector
            else:
                print('Число все еще в запрещенном диапазоне')
                return None
        elif inform_errors >= osn_inform_mult or inform_error_vector_checking > max_error_correction:
            control_error_vector = [(int((A/osn_inform_mult))*osn_inform_mult) % osn[i] for i in range(len(osn))]
            control_error_vector_checking = sum(1 for x in control_error_vector if x != 0)
            if control_error_vector_checking <= max_error_correction:
                print('Ошибка по контрольным основаниям')
                print('Вектор ошибки = ', control_error_vector)
                right_modular_vector = [(chisl[i] - control_error_vector[i]) % osn[i] for i in range(len(osn))]
                print('Правильное число = ', right_modular_vector)
                print('Проверим его на правильность')
                A = check_number(osn, right_modular_vector)
                print('Число в позиционной системе равно = ', A)
                if A < osn_inform_mult:
                    print('Число находится в разрешенном диапазоне, ошибка исправлена правильно')
                    return right_modular_vector
                else:
                    print('Число все еще в запрещенном диапазоне')
                    return None
            elif control_error_vector_checking > max_error_correction and max_error_correction !=1:
                control_modules_for_checking = []
                for combination in combinations(osn_kontrol, max_error_correction):
                    control_osn_combiation_mult = math.prod(combination)
                    control_modules_for_checking.append(control_osn_combiation_mult)
                for combination in control_modules_for_checking:
                    combined_error_vector = [(int((A/combination))*combination) % osn[i] for i in range(len(osn))]
                    combined_error_vector_checking = sum(1 for x in combined_error_vector if x != 0)
                    if combined_error_vector == chisl.tolist() and combined_error_vector_checking > max_error_correction:
                        combined_error_vector = [((int(A/combination)*(combination)) - combination) % osn[i] for i in range(len(osn))]
                        combined_error_vector_checking = sum(1 for x in combined_error_vector if x != 0)
                    elif combined_error_vector == chisl.tolist() and combined_error_vector_checking <= max_error_correction:
                        combined_error_vector = [((int(A/combination)*(combination))) % osn[i] for i in range(len(osn))]
                        combined_error_vector_checking = sum(1 for x in combined_error_vector if x != 0)
                    if combined_error_vector_checking <= max_error_correction and (A - int(A/combination)*combination) < osn_inform_mult:
                        print('Ошибка по информационным и контрольным основаниям')
                        print('Вектор ошибки: ', combined_error_vector)
                        right_modular_vector = [(chisl[i] - combined_error_vector[i]) % osn[i] for i in range(len(osn))]
                        print('Правильное число = ', right_modular_vector)
                        print('Проверим его на правильность')
                        check_number_orth_bas = check_number(osn, right_modular_vector)
                        print('Число в позиционной системе равно = ', check_number)
                        if check_number_orth_bas < osn_inform_mult:
                            print('Число находится в разрешенном диапазоне, ошибка исправлена правильно')
                            return right_modular_vector
                        else:
                            print('Число все еще в запрещенном диапазоне')
                            continue
                    if combined_error_vector_checking > max_error_correction and control_modules_for_checking.index(combination) == k-1:
                        print('Ошибку не удается исправить в заданной модулярной системе')
                        return None
            elif control_error_vector_checking > max_error_correction and max_error_correction !=1:
                print('Ошибку не удается исправить в заданной модулярной системе')
                return None
    return None

def check_with_projections(osn_inform, osn_kontrol, chisl, osn, max_error_correction):
    """Проверяем число методом проекций"""
    osn_kolvo = len(osn)
    work_diapazon = numpy.prod(osn_inform.astype(object))
    k = len(osn_kontrol)
    max_errors_to_correct = max_error_correction
    osn_int = osn.astype(int)
    for num_errors in range(1, max_errors_to_correct + 1):
        for error_positions in combinations(range(osn_kolvo), num_errors):
            proj_osn = numpy.delete(osn_int, error_positions)
            proj_chisl = numpy.delete(chisl, error_positions)
            B_proj = calculate_orthogonal_bases(proj_osn)
            A_proj = numpy.mod(numpy.sum(proj_chisl * B_proj), numpy.prod(proj_osn.astype(object)))
            if A_proj < work_diapazon:
                B_full = calculate_orthogonal_bases(osn_int)
                M = numpy.zeros((num_errors, num_errors), dtype=float)
                for i, pos in enumerate(error_positions):
                    for j, p in enumerate(error_positions):
                        M[i,j] = float(B_full[pos] % osn_int[p])
                b = numpy.zeros(num_errors, dtype=float)
                for i, pos in enumerate(error_positions):
                    s = sum(chisl[p] * B_full[p] for p in range(osn_kolvo) if p not in error_positions)
                    b[i] = float((A_proj - s) % osn_int[pos])
                try:
                    solutions = numpy.linalg.solve(M, b)
                    solutions = numpy.mod(numpy.round(solutions).astype(int), osn_int[list(error_positions)])
                    corrected = chisl.copy()
                    for i, pos in enumerate(error_positions):
                        corrected[pos] = int(solutions[i])
                    A_corr = numpy.mod(numpy.sum(corrected * B_full), numpy.prod(osn_int.astype(object)))
                    if A_corr < work_diapazon:
                        return corrected.tolist()
                except numpy.linalg.LinAlgError:
                    continue
    for num_errors in range(1, k + 1):
        all_valid = True
        for error_positions in combinations(range(osn_kolvo), num_errors):
            temp_chisl = chisl.copy()
            for pos in error_positions:
                temp_chisl[pos] = (temp_chisl[pos] + 1) % osn_int[pos]
            B_full = calculate_orthogonal_bases(osn_int)
            A_err = numpy.mod(numpy.sum(temp_chisl * B_full), numpy.prod(osn_int.astype(object)))
            if A_err < work_diapazon:
                all_valid = False
                break
        if all_valid:
            return None
    
    return None

def main():
    #Ввод данных
    osn_inform_kolvo = int(input('Введите количество информационных оснований: '))
    osn_kontrol_kolvo = int(input('Введите количество контрольных оснований: '))
    osn_inform = numpy.array([int(input(f'Введите {i+1} информационное основание: ')) 
                            for i in range(osn_inform_kolvo)], dtype=object)
    osn_kontrol = numpy.array([int(input(f'Введите {i+1} контрольное основание: ')) 
                            for i in range(osn_kontrol_kolvo)], dtype=object)
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    max_error_correction = int(input('Введите максимальную кратность ошибки: '))
    print('max_error_correction = ', max_error_correction)
    #Ввод числа для проверки
    print('\nВведите число для проверки:')
    chisl = numpy.array([int(input(f'Введите {i+1} число: ')) 
                        for i in range(len(osn))], dtype=object)
    #Проверка новым методом
    new_method_right_vector = new_method_decoding(osn_inform, osn_kontrol, chisl, max_error_correction)
    #Проверка методом проекций
    corrections_proj = check_with_projections(osn_inform, osn_kontrol, chisl, osn, max_error_correction)
    #Вывод результатов
    print("\nРезультаты нового метода")
    print("Исправленное число:", new_method_right_vector)
    print("\nРезультаты метода проекций")
    print("Исправленное число:", corrections_proj)
    # Сравнение результатов
    if corrections_proj and numpy.array_equal(new_method_right_vector, corrections_proj):
        print("\nОба метода дали одинаковый результат!")
    else:
        print("\nРезультаты методов различаются!")
if __name__ == "__main__":
    main()
