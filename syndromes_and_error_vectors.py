import math
import itertools
from collections import defaultdict, Counter
import concurrent.futures
import os
import sys

class ModularSystemAnalyzer:
    def __init__(self, info_bases, control_bases, max_error_multiplicity=None):
        self.info_bases = info_bases
        self.control_bases = control_bases
        self.all_bases = info_bases + control_bases
        self.max_error_multiplicity = max_error_multiplicity
        
        #Вычисляем диапазоны
        self.Pi = math.prod(info_bases)  #Информационный диапазон
        self.Pfull = math.prod(self.all_bases)  #Полный диапазон
        
        #Количество секций
        self.num_sections = self.Pfull // self.Pi
        
        #Генерируем все числа в информационном диапазоне
        self.info_numbers = self._generate_info_numbers()
        
    def _generate_info_numbers(self):
        """Генерирует все числа в информационном диапазоне"""
        numbers = []
        for i in range(self.Pi):
            modular_rep = [i % base for base in self.all_bases]
            numbers.append((i, modular_rep))
        return numbers
    
    def get_modular_representation(self, decimal_num):
        """Возвращает модулярное представление десятичного числа"""
        return [decimal_num % base for base in self.all_bases]
    
    def hamming_distance(self, vec1, vec2):
        """Вычисляет расстояние Хэмминга между двумя векторами"""
        return sum(1 for a, b in zip(vec1, vec2) if a != b)
    
    def get_syndrome(self, section_index):
        """Вычисляет синдром для секции"""
        first_num_in_section = section_index * self.Pi
        modular_rep = self.get_modular_representation(first_num_in_section)
        #Берем последние k элементов (k = количество контрольных модулей)
        return modular_rep[-len(self.control_bases):]
    
    def find_closest_info_number(self, control_number_rep):
        """Находит ближайшее информационное число к заданному контрольному числу"""
        min_distance = float('inf')
        closest_info = None
        
        for info_decimal, info_rep in self.info_numbers:
            distance = self.hamming_distance(control_number_rep, info_rep)
            if distance < min_distance:
                min_distance = distance
                closest_info = (info_decimal, info_rep)
        
        return min_distance, closest_info
    
    def calculate_error_vector(self, control_rep, info_rep):
        """Вычисляет вектор ошибки"""
        error_vector = []
        for c, i, base in zip(control_rep, info_rep, self.all_bases):
            error = (i - c) % base
            error_vector.append(error)
        return error_vector
    
    def get_error_multiplicity(self, error_vector):
        """Вычисляет кратность ошибки (количество ненулевых элементов)"""
        return sum(1 for e in error_vector if e != 0)
    
    def analyze_section(self, section_index):
        """Анализирует одну секцию"""
        start_i = section_index * self.Pi
        end_i = start_i + self.Pi - 1
        
        syndrome = self.get_syndrome(section_index)
        
        section_results = {
            'section_index': section_index,
            'start_i': start_i,
            'end_i': end_i,
            'syndrome': syndrome,
            'numbers': []
        }
        
        error_vectors_in_section = []
        
        #Анализируем все числа в секции
        for i in range(start_i, end_i + 1):
            control_rep = self.get_modular_representation(i)
            
            #Находим ближайшее информационное число
            distance, closest_info = self.find_closest_info_number(control_rep)
            
            if closest_info is not None:
                info_decimal, info_rep = closest_info
                error_vector = self.calculate_error_vector(control_rep, info_rep)
                error_multiplicity = self.get_error_multiplicity(error_vector)
                
                #Проверяем кратность ошибки, если задана максимальная кратность
                if (self.max_error_multiplicity is not None and 
                    error_multiplicity > self.max_error_multiplicity):
                    correction_info = "число не поддается исправлению с заданными параметрами"
                else:
                    correction_info = {
                        'distance': distance,
                        'closest_info_decimal': info_decimal,
                        'closest_info_rep': info_rep,
                        'error_vector': error_vector
                    }
                
                section_results['numbers'].append({
                    'decimal': i,
                    'modular_rep': control_rep,
                    'correction_info': correction_info
                })
                
                if isinstance(correction_info, dict):
                    error_vectors_in_section.append(tuple(error_vector))
        
        #Статистика по ошибочным векторам в секции
        error_vector_stats = Counter(error_vectors_in_section)
        section_results['error_vector_stats'] = error_vector_stats
        section_results['unique_error_vectors'] = len(error_vector_stats)
        
        return section_results
    
    def format_modular_rep(self, rep):
        """Форматирует модулярное представление для вывода"""
        return "[" + ", ".join(map(str, rep)) + "]"
    
    def format_error_info(self, number_data):
        """Форматирует информацию об ошибке для вывода"""
        if number_data['correction_info'] == "число не поддается исправлению с заданными параметрами":
            return "    число не поддается исправлению с заданными параметрами"
        
        info = number_data['correction_info']
        control_rep_str = self.format_modular_rep(number_data['modular_rep'])
        info_rep_str = self.format_modular_rep(info['closest_info_rep'])
        error_vector_str = self.format_modular_rep(info['error_vector'])
        
        lines = [
            f"    Расстояние: {info['distance']}",
            f"    Ближайшее информационное число: i={info['closest_info_decimal']} -> {info_rep_str}",
            f"    Вектор ошибки: {error_vector_str}",
            f"    {control_rep_str} + {error_vector_str} = {info_rep_str}"
        ]
        
        return "\n".join(lines)
    
    def print_section_results(self, section_results, output_file=None):
        """Выводит результаты анализа секции"""
        lines = []
        
        lines.append(f"СЕКЦИЯ {section_results['section_index'] + 1}: i={section_results['start_i']} - i={section_results['end_i']}")
        lines.append(f"Синдром: {section_results['syndrome']}")
        lines.append("")
        
        #Выводим все числа в секции
        for number_data in section_results['numbers']:
            lines.append(f"i={number_data['decimal']}: {self.format_modular_rep(number_data['modular_rep'])}")
            lines.append(self.format_error_info(number_data))
            lines.append("")
        
        #Выводим статистику по секции
        lines.append(f"Количество чисел в секции: {len(section_results['numbers'])}")
        lines.append(f"Уникальных ошибочных векторов: {section_results['unique_error_vectors']}")
        lines.append("=" * 80)
        
        #Записываем в файл или выводим на экран
        result_text = "\n".join(lines)
        if output_file:
            output_file.write(result_text + "\n")
        else:
            print(result_text)
        
        return result_text
    
    def analyze_all_sections_parallel(self, max_workers=None, max_sections=None):
        """Анализирует все секции параллельно"""
        sections_to_analyze = range(self.num_sections)
        if max_sections:
            sections_to_analyze = range(min(max_sections, self.num_sections))
        
        all_results = []
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            #Создаем отдельные экземпляры анализатора для каждого процесса
            future_to_section = {
                executor.submit(self.analyze_section, section_idx): section_idx 
                for section_idx in sections_to_analyze
            }
            
            for future in concurrent.futures.as_completed(future_to_section):
                section_idx = future_to_section[future]
                try:
                    result = future.result()
                    all_results.append(result)
                except Exception as e:
                    print(f"Ошибка при анализе секции {section_idx}: {e}")
        
        #Сортируем результаты по номеру секции
        all_results.sort(key=lambda x: x['section_index'])
        return all_results
    
    def generate_final_statistics(self, all_results):
        """Генерирует итоговую статистику"""
        total_sections = len(all_results)
        total_numbers = sum(len(section['numbers']) for section in all_results)
        
        #Собираем все ошибочные вектора из всех секций
        all_error_vectors = []
        error_vector_details = defaultdict(lambda: {'count': 0, 'sections': set(), 'decimal_ranges': []})
        
        for section in all_results:
            section_idx = section['section_index']
            
            for number_data in section['numbers']:
                if isinstance(number_data['correction_info'], dict):
                    error_vector = tuple(number_data['correction_info']['error_vector'])
                    all_error_vectors.append(error_vector)
                    
                    #Обновляем детальную информацию
                    error_vector_details[error_vector]['count'] += 1
                    error_vector_details[error_vector]['sections'].add(section_idx)
                    error_vector_details[error_vector]['decimal_ranges'].append(number_data['decimal'])
        
        #Сортируем ошибочные вектора по минимальному десятичному числу в диапазоне
        def vector_sort_key(vec):
            details = error_vector_details[vec]
            if details['decimal_ranges']:
                return min(details['decimal_ranges'])
            return float('inf')
        
        sorted_vectors = sorted(set(all_error_vectors), key=vector_sort_key)
        total_unique_vectors = len(sorted_vectors)
        
        #Формируем статистику
        stats_lines = []
        stats_lines.append("=" * 100)
        stats_lines.append("ИТОГОВАЯ СТАТИСТИКА")
        stats_lines.append("=" * 100)
        stats_lines.append(f"Всего проанализировано секций: {total_sections}")
        stats_lines.append(f"Всего проанализировано чисел: {total_numbers}")
        stats_lines.append(f"Всего уникальных ошибочных векторов: {total_unique_vectors}")
        stats_lines.append("")
        stats_lines.append("Ошибочные вектора (отсортировано по диапазону десятичных чисел):")
        
        for i, vec in enumerate(sorted_vectors, 1):
            details = error_vector_details[vec]
            decimal_ranges = details['decimal_ranges']
            min_dec = min(decimal_ranges) if decimal_ranges else 0
            max_dec = max(decimal_ranges) if decimal_ranges else 0
            
            sections_str = ", ".join(map(str, sorted(details['sections'])))
            percentage = (details['count'] / total_numbers * 100) if total_numbers > 0 else 0
            
            #Получаем модулярные представления для начала и конца диапазона
            min_modular_rep = self.get_modular_representation(min_dec)
            max_modular_rep = self.get_modular_representation(max_dec)
            
            stats_lines.append(f"    {i}. {list(vec)}: {details['count']} раз ({percentage:.1f}%)")
            stats_lines.append(f"        Диапазон десятичных чисел: {min_dec} - {max_dec}")
            stats_lines.append(f"        {min_dec} -> {self.format_modular_rep(min_modular_rep)}")
            stats_lines.append(f"        {max_dec} -> {self.format_modular_rep(max_modular_rep)}")
            stats_lines.append(f"        Секции: {sections_str}")
            stats_lines.append("")
        
        return "\n".join(stats_lines)


def main():
    print("Анализатор модулярной системы счисления")
    print("=" * 50)
    
    #Ввод параметров от пользователя
    try:
        num_info = int(input("Введите количество информационных модулей: "))
        num_control = int(input("Введите количество контрольных модулей: "))
        
        print(f"Введите {num_info} информационных модулей:")
        info_bases = []
        for i in range(num_info):
            base = int(input(f"Информационный модуль {i+1}: "))
            info_bases.append(base)
        
        print(f"Введите {num_control} контрольных модулей:")
        control_bases = []
        for i in range(num_control):
            base = int(input(f"Контрольный модуль {i+1}: "))
            control_bases.append(base)
        
        #Опциональный параметр максимальной кратности ошибки
        max_error_input = input("Введите максимальную кратность ошибки (или нажмите Enter для пропуска): ")
        max_error_multiplicity = int(max_error_input) if max_error_input.strip() else None
        
        max_sections_input = input("Введите максимальное количество секций для анализа (или нажмите Enter для всех): ")
        max_sections = int(max_sections_input) if max_sections_input.strip() else None
        
        workers_input = input("Введите количество рабочих процессов (или нажмите Enter для автоматического определения): ")
        max_workers = int(workers_input) if workers_input.strip() else os.cpu_count()
        
    except ValueError as e:
        print(f"Ошибка ввода: {e}")
        return
    
    #Создаем анализатор
    analyzer = ModularSystemAnalyzer(info_bases, control_bases, max_error_multiplicity)
    
    print("\nПараметры системы:")
    print(f"Информационные модули: {info_bases}")
    print(f"Контрольные модули: {control_bases}")
    print(f"Информационный диапазон (Pi): {analyzer.Pi}")
    print(f"Полный диапазон (Pfull): {analyzer.Pfull}")
    print(f"Количество секций: {analyzer.num_sections}")
    print(f"Максимальная кратность ошибки: {max_error_multiplicity}")
    
    #Предупреждение о размере файла
    estimated_size_mb = (analyzer.Pfull * 0.1) / 1024 / 1024  # Примерная оценка
    print(f"\nПримерный размер выходного файла: {estimated_size_mb:.2f} MB")
    
    file_output = input("Записать результаты в файл? (y/n): ").lower().strip() == 'y'
    screen_output = input("Вывести результаты на экран? (y/n): ").lower().strip() == 'y'
    
    output_file = None
    if file_output:
        filename = input("Введите имя файла для сохранения: ")
        try:
            output_file = open(filename, 'w', encoding='utf-8')
            print(f"Результаты будут сохранены в файл: {filename}")
        except Exception as e:
            print(f"Ошибка при создании файла: {e}")
            file_output = False
    
    #Анализ секций
    print("\nНачинаем анализ секций...")
    all_results = analyzer.analyze_all_sections_parallel(
        max_workers=max_workers, 
        max_sections=max_sections
    )
    
    #Вывод результатов
    print("\nРезультаты анализа:")
    for section_results in all_results:
        analyzer.print_section_results(section_results, output_file if file_output else None)
    
    #Итоговая статистика
    final_stats = analyzer.generate_final_statistics(all_results)
    
    if output_file:
        output_file.write(final_stats + "\n")
        output_file.close()
        print(f"\nРезультаты сохранены в файл: {filename}")
    
    if screen_output:
        print(final_stats)
    
    print("\nАнализ завершен!")


if __name__ == "__main__":
    main()
