import os
import subprocess

# Caminho para o diretório de instâncias
instances_dir = "Instances/csarp/"

# Caminho para o executável
bin_executable = "./bin/bin.out"

# Argumentos fixos
arg_1 = "1AD"
bundle_output_file = "file1.txt"
node_output_file = "file2.txt"

# Pasta de destino para os resultados comparados
compare_solutions_dir = "compareResults/"

# Cria a pasta se ela não existir
os.makedirs(compare_solutions_dir, exist_ok=True)

# Função para extrair tempo e LB
def extract_time_and_lb(file_path):
    time_value = None
    lb_value = None
    with open(file_path, "r") as file:
        for line in file:
            if "Total Time:" in line:
                time_value = float(line.split("Total Time:")[1].strip())
            if "LB:" in line:
                lb_value = float(line.split("LB:")[1].strip())
    return time_value, lb_value

# Função para executar o comando e redirecionar a saída para um arquivo
def run_command(instance, arg_2, output_file):
    instance_path = os.path.join(instances_dir, instance)
    command = f"{bin_executable} {instance_path} {arg_1} {arg_2} >> {output_file}"
    subprocess.run(command, shell=True)

# Função principal para rodar a instância em ambos os casos
for instance in os.listdir(instances_dir):
    if instance.endswith(".txt"):  # Verifica se é um arquivo .txt
        # Executa com "bundle2" e armazena a saída em file1.txt
        run_command(instance, "bundle2", bundle_output_file)

        # Executa com "node" e armazena a saída em file2.txt
        run_command(instance, "node", node_output_file)

        # Extrai resultados do file1.txt (bundle2) e file2.txt (node)
        bundle_time, bundle_lb = extract_time_and_lb(bundle_output_file)
        node_time, node_lb = extract_time_and_lb(node_output_file)

        # Armazena os resultados combinados na pasta compareResults/
        compare_result_path = os.path.join(compare_solutions_dir, f"{instance}_compare.txt")

        with open(compare_result_path, "w") as compare_file:
            compare_file.write(f"Instância: {instance}\n")
            compare_file.write("==== Bundle2 ====\n")
            compare_file.write(f"LB: {bundle_lb}\n")
            compare_file.write(f"Total Time: {bundle_time} segundos\n\n")

            compare_file.write("==== Node ====\n")
            compare_file.write(f"LB: {node_lb}\n")
            compare_file.write(f"Total Time: {node_time} segundos\n\n")

            # Calcula a diferença de tempo entre bundle2 e node
            time_difference = bundle_time - node_time
            compare_file.write(f"==== Comparação ====\n")
            compare_file.write(f"Diferença de tempo (Bundle2 - Node): {time_difference} segundos\n")

        # Exibe os resultados no terminal para cada instância
        print(f"Instância: {instance}")
        print(f"Bundle2 - LB: {bundle_lb}, Tempo: {bundle_time}s")
        print(f"Node - LB: {node_lb}, Tempo: {node_time}s")
        print(f"Diferença de tempo: {time_difference}s\n")

print("Execução finalizada e resultados comparados em compareResults/!")
