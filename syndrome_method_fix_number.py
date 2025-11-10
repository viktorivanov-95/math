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

def solve_linear_equation(equation_str, var_name, known_vars=None, mod=None):
    if known_vars is None:
        known_vars = {}
    
    for var, value in known_vars.items():
        equation_str = equation_str.replace(var, str(value))
    
    if mod is not None:
        for x in chain([0], (x*sign for x in count(1) for sign in (1, -1))):
            x_mod = x % mod
            expr = equation_str.replace(var_name, str(x_mod))
            try:
                value = eval(expr, {'__builtins__': None}, known_vars)
                if value % mod == 0:
                    return x_mod
            except:
                continue
    return None

#Вычисление синдромов для обнаруженных ошибок
def syndrome_decoding(osn_inform, osn_kontrol, chisl):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    k = len(osn_kontrol)
    n = len(osn_inform)
    n_plus_k = len(osn)
    B = calculate_orthogonal_bases(osn)
    B = numpy.array(B, dtype=object)
    print('B = ', B)    
    B_modular = numpy.zeros((n_plus_k, n_plus_k), dtype=object)
    for i in range(n_plus_k):
        B_modular[i][:] = [b % int(osn[i]) for b in B]
    print('B_modular = ')
    print(B_modular)
    c = numpy.zeros((n_plus_k-1, n_plus_k), dtype=object)
    for i in range(n_plus_k-1):
        for j in range(n_plus_k-1):
            try:
                c[j][i] = pow(int(osn[i]), -1, int(osn[j+1]))
            except ValueError:
                continue
    
    for i in range(n_plus_k-2):
        c[i][i+1:] = 0
    
    B_opss = numpy.zeros((n_plus_k, n_plus_k), dtype=object)
    for k in range(n_plus_k):
        B_opss[k][0] = B_modular[k][0] % osn[0]
        if n_plus_k > 1:
            B_opss[k][1] = ((B_modular[k][1] - B_opss[k][0]) * c[0][0]) % osn[1]
        for i in range(2, n_plus_k):
            temp = (B_modular[k][i] - B_opss[k][0]) * c[i-1][0]
            for j in range(1, i):
                temp = (temp - B_opss[k][j]) * c[i-1][j]
            B_opss[k][i] = temp % osn[i]
    print('B_opss = ')
    print(B_opss)
    B_opss_chisl = numpy.zeros((n_plus_k, n_plus_k), dtype=object)
    for i in range(n_plus_k):
        for j in range(n_plus_k):
            B_opss_chisl[i][j] = (B_opss[i][j] * chisl[i]) if i < n else B_opss[i][j]
    print('B_opss_chisl = ')
    print(B_opss_chisl)
    opss = numpy.zeros(n_plus_k, dtype=object)
    for j in range(n_plus_k):
        for i in range(n_plus_k):
            if i < n:
                opss[j] += B_opss_chisl[i][j]
            else:
                if B_opss_chisl[i][j] != 0:
                    var_name = f"x{i-n+1}"
                    term = f"{B_opss_chisl[i][j]}*{var_name}"
                    opss[j] = f"{opss[j]} + {term}" if opss[j] != 0 else term
                    
    opss_total = numpy.zeros(n_plus_k, dtype=object)
    perenos = 0
    
    for i in range(n_plus_k):
        if isinstance(opss[i], (int, numpy.integer)):
            total = int(opss[i]) + int(perenos)
            opss_total[i] = total % osn[i]
            perenos = total // osn[i]
        else:
            opss_total[i] = f"({opss[i]} + {perenos})"
            perenos = f"({opss[i]} + {perenos}) // {osn[i]}"

    print('opss_total = ')
    print(opss_total)

    syndromes = []
    known_vars = {}
    
    for i in range(n, n_plus_k):
        mod = int(osn[i])
        equation = str(opss_total[i])
        print('equation = ')
        print(equation, ' = 0')
        var_name = f"x{i-n+1}"
        
        solution = solve_linear_equation(equation, var_name, known_vars, mod)

        print('solution = ', solution)
        
        if solution is not None:
            known_vars[var_name] = solution
            syndrome = (int(chisl[i]) - solution) % mod
            print('syndrome = (', int(chisl[i]), ' - ', solution, ') % ', mod)
            print('syndrome = ', syndrome)
            syndromes.append(syndrome)
        else:
            syndromes.append(None)

    return syndromes

#Исправление ошибок на основе синдромов
def locate_and_correct_errors(osn_inform, osn_kontrol, chisl, syndromes):
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    k = len(osn_kontrol)
    n = len(osn_inform)
    work_diapazon = math.prod([int(x) for x in osn_inform])
    print('syndromes = ', syndromes)
    all_zero = all(s == 0 or s is None for s in syndromes)
    all_non_zero = all(s != 0 and s is not None for s in syndromes)
    some_non_zero = any(s != 0 and s is not None for s in syndromes) and not all_non_zero
    
    if all_zero:
        print("\nВсе синдромы равны нулю - ошибок нет")
        return chisl
    
    elif all_non_zero:
        print("\nВсе синдромы ненулевые - ошибка в информационных основаниях")
        P = work_diapazon
        P_i = [P // int(x) for x in osn_inform]
        
        inv_P_i = []
        for i in range(n):
            inv_row = []
            for j in range(k):
                try:
                    inv = pow(int(P_i[i]), -1, int(osn[n+j]))
                    inv_row.append(inv)
                except ValueError:
                    inv_row.append(None)
            inv_P_i.append(inv_row)

        lambda_values_minus = []
        lambda_values_plus = []
        for i in range(n):
            lambda_row_minus = []
            lambda_row_plus = []
            for j in range(k):
                if inv_P_i[i][j] is not None and syndromes[j] is not None:
                    lambda_val_minus = (-syndromes[j] * inv_P_i[i][j]) % osn[n+j]
                    print('lambda_val_minus = (', -syndromes[j], '*', inv_P_i[i][j], ')%', osn[n+j], ' = ', lambda_val_minus)
                    lambda_val_plus = (syndromes[j] * inv_P_i[i][j]) % osn[n+j]
                    print('lambda_val_plus = (', syndromes[j], '*', inv_P_i[i][j], ')%', osn[n+j], ' = ', lambda_val_plus)
                    lambda_row_minus.append(lambda_val_minus)
                    lambda_row_plus.append(lambda_val_plus)
                else:
                    lambda_row_minus.append(None)
                    lambda_row_plus.append(None)
            lambda_values_minus.append(lambda_row_minus)
            lambda_values_plus.append(lambda_row_plus)

        print('lambda_values_minus = ', lambda_values_minus)
        print('lambda_values_plus = ', lambda_values_plus)

        error_pos_minus = None
        error_value_minus = None
        error_pos_plus = None
        error_value_plus = None
        
        for i in range(n):
            unique_lambdas_minus = set(val for val in lambda_values_minus[i] if val is not None)
            if len(unique_lambdas_minus) == 1 and unique_lambdas_minus:
                error_pos_minus = i
                error_value_minus = (unique_lambdas_minus.pop() * P_i[i]) % osn_inform[i]
                break
            
        for i in range(n):
            unique_lambdas_plus = set(val for val in lambda_values_plus[i] if val is not None)
            if len(unique_lambdas_plus) == 1 and unique_lambdas_plus:
                error_pos_plus = i
                error_value_plus = (unique_lambdas_plus.pop() * P_i[i]) % osn_inform[i]
                break

        corrected_minus = None
        corrected_plus = None

        if error_pos_minus is not None:
            corrected_minus = numpy.array(chisl, dtype=object)
            corrected_minus[error_pos_minus] = (corrected_minus[error_pos_minus] - error_value_minus) % osn_inform[error_pos_minus]
            print('corrected_minus = ', corrected_minus)
        
        if error_pos_plus is not None:
            corrected_plus = numpy.array(chisl, dtype=object)
            corrected_plus[error_pos_plus] = (corrected_plus[error_pos_plus] + error_value_plus) % osn_inform[error_pos_plus]
            print('corrected_plus = ', corrected_plus)
            
        valid_corrections = []
        
        if corrected_minus is not None:
            A_minus = check_number(osn, corrected_minus)
            if A_minus < work_diapazon:
                valid_corrections.append(('minus', corrected_minus, error_value_minus, error_pos_minus))
        
        if corrected_plus is not None:
            A_plus = check_number(osn, corrected_plus)
            if A_plus < work_diapazon:
                valid_corrections.append(('plus', corrected_plus, error_value_plus, error_pos_plus))
        
        if valid_corrections:
            correction_type, corrected, error_value, error_pos = valid_corrections[0]
            print(f"Исправлено {'вычитанием' if correction_type == 'minus' else 'прибавлением'} {error_value} {'из' if correction_type == 'minus' else 'к'} позиции {error_pos}")
            return corrected
        else:
            print("Не удалось исправить ошибку")
            return chisl
    
    elif some_non_zero:
        print("\nОбнаружена ошибка в контрольных основаниях")
        corrected = numpy.array(chisl, dtype=object)
        for i in range(k):
            if syndromes[i] != 0 and syndromes[i] is not None:
                B = calculate_orthogonal_bases(osn)
                A = sum(int(chisl[j]) * int(B[j]) for j in range(n)) % work_diapazon
                corrected[n+i] = A % osn[n+i]
                print(f"Исправлено контрольное основание {osn[n+i]}: было {chisl[n+i]}, стало {corrected[n+i]}")
        return corrected
    
    else:
        print("\nНе удалось определить тип ошибки")
        return chisl

def check_with_projections(osn_inform, osn_kontrol, chisl, osn):
    """Проверяет число методом проекций."""
    osn_kolvo = len(osn)
    work_diapazon = numpy.prod(osn_inform.astype(object))
    k = len(osn_kontrol)
    max_errors_to_correct = ceil(k / 2)
    
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
                        corrected[pos] = solutions[i]
                    
                    A_corr = numpy.mod(numpy.sum(corrected * B_full), numpy.prod(osn_int.astype(object)))
                    if A_corr < work_diapazon:
                        corrections = [(pos, int(solutions[i])) for i, pos in enumerate(error_positions)]
                        return corrections, f"Исправлено {num_errors} ошибок"
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
            return None, f"Обнаружено {num_errors} ошибок (исправление невозможно)"
    
    return None, "Не удалось обнаружить ошибки (возможно, их больше чем k)"

def main():
    #Ввод данных
    osn_inform_kolvo = int(input('Введите количество информационных оснований: '))
    osn_kontrol_kolvo = int(input('Введите количество контрольных оснований: '))
    osn_inform = numpy.array([int(input(f'Введите {i+1} информационное основание: ')) 
                            for i in range(osn_inform_kolvo)], dtype=object)
    osn_kontrol = numpy.array([int(input(f'Введите {i+1} контрольное основание: ')) 
                            for i in range(osn_kontrol_kolvo)], dtype=object)
    osn = numpy.concatenate((osn_inform, osn_kontrol))
    
    #Ввод числа для проверки
    print('\nВведите число для проверки:')
    chisl = numpy.array([int(input(f'Введите {i+1} число: ')) 
                        for i in range(len(osn))], dtype=object)
    
    #Проверка методом синдромного декодирования
    syndromes = syndrome_decoding(osn_inform, osn_kontrol, chisl)
    corrected_syndrom = locate_and_correct_errors(osn_inform, osn_kontrol, chisl, syndromes)
    
    #Проверка методом проекций
    corrections_proj, message_proj = check_with_projections(osn_inform, osn_kontrol, chisl, osn)
    
    #Вывод результатов
    print("\nРезультаты синдромного декодирования")
    print("Исправленное число:", corrected_syndrom)
    
    print("\nРезультаты метода проекций")
    if corrections_proj:
        corrected_proj = chisl.copy()
        for pos, correct_value in corrections_proj:
            corrected_proj[pos] = correct_value
        print("Исправленное число:", corrected_proj)
    else:
        print(message_proj)
    
    #Сравнение результатов
    if corrections_proj and numpy.array_equal(corrected_syndrom, corrected_proj):
        print("\nОба метода дали одинаковый результат!")
    else:
        print("\nРезультаты методов различаются!")

if __name__ == "__main__":
    main()
