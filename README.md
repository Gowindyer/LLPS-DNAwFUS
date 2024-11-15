LLPS_FUSandDNA是基于openMM 开发来模拟fus蛋白和DNA液液相分离的toy模型。fus蛋白和DNA皆由一条单链表示，长短可以根据需求自行设定。主要文件介绍如下：

# Force field
里面是力场文件，force.xml,主要包括原子类型、残基类型和力。力主要是Harmonica bond force以及两个自定义非键力排斥力, fus蛋白链间的L—J力,fus蛋白和DNA之间的L-J，详细见源文件。

# gene_system_pdb
该文件夹是用来创建两个小分子(DNA and protein fus)和混合系统的,会生成pdb和psf文件。构建系统的参数文件simulation_parameters.json

生成系统命令：./run_construct_system.sh

# run_simulation.py
执行模拟的文件，模拟参数文件可以通过修改simulaiton_parameters.json的指定内容来修改。注意，当前脚本运行在CPU上,如需GPU，请修改simulation_parameters.json里的计算平台选项为CUDA。
执行命令：python run_simulation.py
CPU指定核数(如1核)： OPENMM_CPU_THREADS=1 python run_simulation.py

# 软件包
openmm: 生成两种小分子，运行分子模拟
Packmol: 构建混合系统.
json: 加载参数文件.
等

2024.4.30
* 增加nonbonded_terms 函数
* 力场创建方式由基于xml构建改为基于python脚本函数构建。
* 力场参数(排斥力和lj势)改为由从params里的ex_epsilon.csv、ex_sigma.csv、lj_epsilon.csv和lj_sigma.csv获得。原json文件不再保有力场相关参数。

2024.1.26
更正generate system pdb部分代码：主要包括有
* 修正run_construct_system.sh的路径名
* 修正DNA，protein_fus中construct_polymer.py中的generate PSF 文件的代码
* 修正pack mixed system的gene_psf.py的原子名错误：更正：'NATOME'为'NATOM'
