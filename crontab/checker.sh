source ~/.bashrc
conda activate rmg3
python /home/hpc/ypli/jianyi/test_pygsm/AutomaticReactionDiscovery/database/checker.py
conda deactivate
echo $(date +%Y-%m-%d:%H:%M:%S)