import os
import subprocess

# Caminho para o diretório de instâncias
instances_dir = "Instances/csarp/"

# Caminho para o executável
bin_executable = "./bin/bin.out"

# Argumentos fixos
arg_1 = "1AD"
arg_2 = "node"
output_file = "nodeSolutions.txt"

# Percorre todos os arquivos no diretório de instâncias
for instance in os.listdir(instances_dir):
    if instance.endswith(".txt"):  # Verifica se o arquivo é um .txt
        instance_path = os.path.join(instances_dir, instance)  # Cria o caminho completo
        
        # Comando a ser executado
        command = f"{bin_executable} {instance_path} {arg_1} {arg_2} >> {output_file}"
        
        # Executa o comando no terminal
        subprocess.run(command, shell=True)

print("Comandos executados para todos os arquivos .txt!")
