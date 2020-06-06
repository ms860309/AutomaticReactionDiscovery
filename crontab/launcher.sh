export BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../AutomaticReactionDiscovery" && pwd )"
export PYTHONPATH=$PYTHONPATH:$BASE_DIR
export PATH=$HOME/anaconda3/bin:$PATH

source ~/.bashrc
conda activate rmg3
echo $(date +%Y-%m-%d:%H:%M:%S)
python $BASE_DIR/database/launcher.py
python $BASE_DIR/database/network_launcher.py
conda deactivate