import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import gc
import signal
import psutil

# 定义主目录和相关参数
root_dir = '/data1/bioinfo/zhoutong/data/RASP_tmp/50_PRJNA608297/plus_minus/link/HEK293N_Homo_sapiens'
output_dir = '/data1/bioinfo/zhoutong/data/RASP_tmp/50_PRJNA608297/plus_minus/genome_model_input'
command = '/home/bioinfo/06_dev/Genome_Wide_merge/bin/stone/nodepengsearch/target/release/nodepengsearch'
length_file = '/data1/bioinfo/YSZ/fasta2len/plus_minus_0929.len'
threads = 100

# 定义执行命令的函数
def run_nodepengsearch(input_file):
    output_file = os.path.join(output_dir, os.path.basename(input_file).replace('.txt', '_test.csv'))

    cmd = [
        command, '--len', length_file, '--strand', '+', '--merged', input_file,
        '--output', output_file, '--threads', str(threads)
    ]

    try:
        # 执行命令并等待完成
        subprocess.run(cmd, check=True)
        print(f'处理完成: {input_file}')
    except subprocess.CalledProcessError as e:
        print(f'处理 {input_file} 时出错: {e}')
    finally:
        # 强制进行垃圾回收以释放资源    
        gc.collect()

# 捕获中断信号并释放资源
def signal_handler(sig, frame):
    print("中断信号接收，正在终止所有任务...")
    # 在接收中断时清理资源
    gc.collect()
    os._exit(1)

# 注册信号处理器以捕捉中断信号 (如 Ctrl+C)
signal.signal(signal.SIGINT, signal_handler)

# 遍历主目录下的所有子文件夹，查找以 merged.txt 结尾的文件
def find_merged_files(root_dir):
    merged_files = []
    for dirpath, _, filenames in os.walk(root_dir):
        for file in filenames:
            if file.endswith('merged.txt'):
                merged_files.append(os.path.join(dirpath, file))
    return merged_files

# 获取所有以 merged.txt 结尾的文件的路径
merged_files = find_merged_files(root_dir)

# 定义最大并行任务数
max_parallel_tasks = 15

# 使用线程池并行执行命令，捕捉任务完成或中断
with ThreadPoolExecutor(max_workers=max_parallel_tasks) as executor:
    # 提交任务
    future_to_file = {executor.submit(run_nodepengsearch, file): file for file in merged_files}
    
    try:
        # 逐个等待任务完成，并确保每个任务的内存被释放
        for future in as_completed(future_to_file):
            input_file = future_to_file[future]
            try:
                future.result()
            except Exception as exc:
                print(f'{input_file} 处理时发生异常: {exc}')
            finally:
                # 执行完成后强制进行垃圾回收
                gc.collect()

    except KeyboardInterrupt:
        print("用户中断，正在清理并释放资源...")
        # 在 KeyboardInterrupt 时清理线程池并退出
        executor.shutdown(wait=False)
        gc.collect()
        os._exit(1)

print("所有文件处理完成。")

# 监控内存占用
def monitor_memory():
    process = psutil.Process(os.getpid())
    print(f"当前内存使用: {process.memory_info().rss / (1024 ** 2)} MB")
