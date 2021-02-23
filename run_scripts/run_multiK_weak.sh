#nohup bash -c "matlab -nodesktop < models/multiK_weak_enforcement_n0_667.m > results/multiK_weak_enforcement_n0_667.out; \
#matlab -nodesktop < models/multiK_weak_enforcement_linearized_cov_n0_667.m > results/multiK_weak_enforcement_linearized_cov_n0_667.out; \
#matlab -nodesktop < models/multiK_weak_enforcement_linearized_output_n0_667.m > results/multiK_weak_enforcement_linearized_output_n0_667.out;" > multiK_weak_n0_667.out &

#nohup bash -c "matlab -nodesktop < models/multiK_weak_enforcement_n1_0.m > results/multiK_weak_enforcement_n1_0.out; \
#matlab -nodesktop < models/multiK_weak_enforcement_linearized_cov_n1_0.m > results/multiK_weak_enforcement_linearized_cov_n1_0.out; \
#matlab -nodesktop < models/multiK_weak_enforcement_linearized_output_n1_0.m > results/multiK_weak_enforcement_linearized_output_n1_0.out;" > multiK_weak_n1.out &

#nohup bash -c "matlab -nodesktop < models/multiK_weak_enforcement_n2_5.m > results/multiK_weak_enforcement_n2_5.out; \
#matlab -nodesktop < models/multiK_weak_enforcement_linearized_cov_n2_5.m > results/multiK_weak_enforcement_linearized_cov_n2_5.out; \
#matlab -nodesktop < models/multiK_weak_enforcement_linearized_output_n2_5.m > results/multiK_weak_enforcement_linearized_output_n2_5.out; " > multiK_weak_n2_5.out &

nohup bash -c "matlab -nodesktop < models/multiK_weak_enforcement_n0_667.m > results/multiK_weak_enforcement_n0_667.out;" > multiK_weak_n0_667.out &

nohup bash -c "matlab -nodesktop < models/multiK_weak_enforcement_n2_5.m > results/multiK_weak_enforcement_n2_5.out;" > multiK_weak_n2_5.out &
