import os
import subprocess
import concurrent.futures
import csv


def execute_binary(binary_path, file_path, arg1, arg2):

    index = file_path.rfind("/")
    name = file_path[(index+1):]

    with open(f"bundle4-2MM/{name}.txt", "a") as outfile:
        arguments = [binary_path, file_path, arg1, arg2]

        subprocess.run(arguments, stdout=outfile, stderr=outfile)


def get_instances(file_path):
    instancias_csv  = []
    starts_csv      = []
    # counter = 0

    with open(file_path, mode='r', newline='', encoding='utf-8') as file:
        leitor_csv = csv.reader(file)

        for linha in leitor_csv:
            # print(linha)
            if len(linha) == 2:
                instancias_csv.append(linha[0])
                starts_csv.append(linha[1])
            
            # counter += 1
            # if counter == 5:
            #     break
    
    # for linha in dados:
    #     print(linha)

    return instancias_csv, starts_csv



def run(bin_path, instances, starts, max_threads):
    arg1 = "2MM"
    arg2 = "bundle4"

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_threads) as executor:

        future_to_command = {executor.submit(execute_binary, bin_path, instance, arg1, arg2): instance for instance in instances}

        for future in concurrent.futures.as_completed(future_to_command):
            command = future_to_command[future]
            try:
                future.result()
            except Exception as e:
                print(f"Subprocess for command 'teste' generated an exception: {e}")


tempInst = []
tempStart = []

instsC1, startsC1 = get_instances("arquivos_relacionados_C1.csv")

instsC2, startsC2 = get_instances("arquivos_relacionados_C2.csv")

#instances = [
#    "Instances_M2/csarp/sarp-30-30-C-2.txt",
#    "Instances_M2/sf_data/sfsarp-30-20-1.txt",
#    "Instances_M2/sf_data/sfsarp-30-20-2.txt",
#    "Instances_M2/sf_data/sfsarp-30-30-2.txt",
#    "Instances_M2/sf_data/sfsarp-30-30-3.txt",
#    "Instances_M2/sf_data/sfsarp-30-30-4.txt"
#]
#starts = []

#instances = ["Instances_M2/csarp/sarp-30-30-C-2.txt"]
#starts = []

#run("./bin/bin.out", instances, starts, 1)


#######################################################################
BIN_PATH = "./bin/bin.out"
MAX_THREADS_C1 = 1
MAX_THREADS_C2 = 1

run(BIN_PATH, instsC1, startsC1, MAX_THREADS_C1)
run(BIN_PATH, instsC2, startsC2, MAX_THREADS_C2)
#######################################################################

#print(len(instances))
#for i in instances:
#    print(i)
#print(len(starts))
#for i in starts:
#    print(i)

# for i in range(len(instances)):
#     instances[i] += ".txt"

# instances, starts = get_instances("inst_&_bounds_C2.csv")
# for i in range(len(instances)):
#     instances[i] += ".txt"

# for instance, start in zip(instances, starts):
#     print(f"{instance}, {start}")

# print("instances + bounds:")
# for instance, bound in zip(instances, bounds):
#     print(instance + ", " + str(bound))


######################### Bundle Based #########################
#BUNDLE_BASED_PATH = "./bundle_based_BnC/bin/bin.out"
#BUNDLE_BASED_MAX_THREADS = 1
#instances_bundle = instances.copy()
#for i in range(len(instances_bundle)):
#    instances_bundle[i] = "bundle_based_BnC/" + instances_bundle[i]


#run(BUNDLE_BASED_PATH, "0", False, instances_bundle, starts, BUNDLE_BASED_MAX_THREADS)
# run(BUNDLE_BASED_PATH, "0", True, instances_bundle, starts, BUNDLE_BASED_MAX_THREADS)

#run(BUNDLE_BASED_PATH, "1", False, instances_bundle, starts, BUNDLE_BASED_MAX_THREADS)
#run(BUNDLE_BASED_PATH, "1", True, instances_bundle, starts, BUNDLE_BASED_MAX_THREADS)
################################################################


######################### Node Based #########################
#NODE_BASED_PATH = "./node_based_BnC/bin/"
#NODE_BASED_MAX_THREADS = 20
#instances_node = instances.copy()
#for i in range(len(instances_node)):
#    instances_node[i] = "node_based_BnC/" + instances_node[i]

# LAZY_ROUTE
# LAZY_ROUTE_PATH = NODE_BASED_PATH + "/lazy_route.out"

# run(LAZY_ROUTE_PATH, "0", False, instances_node, bounds, NODE_BASED_MAX_THREADS)
# run(LAZY_ROUTE_PATH, "0", True, instances_node, bounds, NODE_BASED_MAX_THREADS)

# run(LAZY_ROUTE_PATH, "1", False, instances_node, bounds, NODE_BASED_MAX_THREADS)
# run(LAZY_ROUTE_PATH, "1", True, instances_node, bounds, NODE_BASED_MAX_THREADS)


# LAZY_PERMUTATIONS
# LAZY_PERMUTATIONS_PATH = NODE_BASED_PATH + "/lazy_permutations.out"

# run(LAZY_PERMUTATIONS_PATH, "0", False, instances_node, bounds, NODE_BASED_MAX_THREADS)
# run(LAZY_PERMUTATIONS_PATH, "0", True, instances_node, bounds, NODE_BASED_MAX_THREADS)

# run(LAZY_PERMUTATIONS_PATH, "1", False, instances_node, bounds, NODE_BASED_MAX_THREADS)
# run(LAZY_PERMUTATIONS_PATH, "1", True, instances_node, bounds, NODE_BASED_MAX_THREADS)


# TW_BASED
# TW_BASED_PATH = NODE_BASED_PATH + "/tw_with_restrict.out"

# run(TW_BASED_PATH, "0", False, instances_node, bounds, NODE_BASED_MAX_THREADS)
# run(TW_BASED_PATH, "0", True, instances_node, bounds, NODE_BASED_MAX_THREADS)

# run(TW_BASED_PATH, "1", False, instances_node, bounds, NODE_BASED_MAX_THREADS)
# run(TW_BASED_PATH, "1", True, instances_node, bounds, NODE_BASED_MAX_THREADS)


# ASSIGNMENT_BASED
#ASSIGNMENT_BASED = NODE_BASED_PATH + "/assignment_based.out"

#run(ASSIGNMENT_BASED, "0", False, instances_node, starts, NODE_BASED_MAX_THREADS)
#run(ASSIGNMENT_BASED, "0", True, instances_node, starts, NODE_BASED_MAX_THREADS)

#run(ASSIGNMENT_BASED, "1", False, instances_node, starts, NODE_BASED_MAX_THREADS)
#run(ASSIGNMENT_BASED, "1", True, instances_node, starts, NODE_BASED_MAX_THREADS)
##############################################################


######################### Position Based #########################
#POSITION_BASED_PATH = "./position_based_BnC/bin/bin.out"
#POSITION_BASED_MAX_THREADS = 2
#instances_position = instances.copy()
#for i in range(len(instances_position)):
#    instances_position[i] = "position_based_BnC/" + instances_position[i]

#N = 150
#instances_zero = instances_position[N:]
#print(instances_zero[0])
#run(POSITION_BASED_PATH, "0", False, instances_zero, starts, POSITION_BASED_MAX_THREADS)
#run(POSITION_BASED_PATH, "0", True, instances_position, starts, POSITION_BASED_MAX_THREADS)

#run(POSITION_BASED_PATH, "1", False, instances_position, starts, POSITION_BASED_MAX_THREADS)
#run(POSITION_BASED_PATH, "1", True, instances_position, starts, POSITION_BASED_MAX_THREADS)
##################################################################
