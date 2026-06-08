# step 1 
python inference_loop.py #  分别跑batch1,batch2，batch3,....

# step 2
python process_design_results.py  # 处理step1 得到的结果

# step 3
python combined_final_results.py  # 汇总 step2 处理的各个batch的指标文件

# step 4
bash run_rosetta_InterfaceAnalyzer.sh /path/to/final_data   #  执行rosetta InterfaceAnalyzer,输入是step3 汇总得到的final_data目录（目录下的pdb文件）

# step 5
python merge_rosetta_result.py  #  汇总step3 和step4 的指标文件
