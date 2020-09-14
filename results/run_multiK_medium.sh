nohup bash -c "matlab -nodesktop < multiK_medium_enforcement_n0_667.m > multiK_medium_enforcement_n0_667.out; \
matlab -nodesktop < multiK_medium_enforcement_linearized_cov_n0_667.m > multiK_medium_enforcement_linearized_cov_n0_667.out; \
matlab -nodesktop < multiK_medium_enforcement_linearized_output_n0_667.m > multiK_medium_enforcement_linearized_output_n0_667.out;" &

nohup bash -c "matlab -nodesktop < multiK_medium_enforcement_n1_0.m > multiK_medium_enforcement_n1_0.out; \
matlab -nodesktop < multiK_medium_enforcement_linearized_cov_n1_0.m > multiK_medium_enforcement_linearized_cov_n1_0.out; \
matlab -nodesktop < multiK_medium_enforcement_linearized_output_n1_0.m > multiK_medium_enforcement_linearized_output_n1_0.out;" &

nohup bash -c "matlab -nodesktop < multiK_medium_enforcement_n2_5.m > multiK_medium_enforcement_n2_5.out; \
matlab -nodesktop < multiK_medium_enforcement_linearized_cov_n2_5.m > multiK_medium_enforcement_linearized_cov_n2_5.out; \
matlab -nodesktop < multiK_medium_enforcement_linearized_output_n2_5.m > multiK_medium_enforcement_linearized_output_n2_5.out; " &
