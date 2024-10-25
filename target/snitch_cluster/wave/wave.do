onerror {resume}
quietly WaveActivateNextPane {} 0
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/clk_i}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/rst_i}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/trace_port_o}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/inp_qaddr_i}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/inp_qid_i}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/inp_qdata_op_i}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/inp_qdata_arga_i}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/inp_qdata_argb_i}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/inp_qdata_argc_i}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/inp_qvalid_i}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/inp_qready_o}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/oup_qaddr_o}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/oup_qid_o}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/oup_qdata_arga_o}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/oup_qdata_argb_o}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/oup_qdata_argc_o}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/oup_qdata_repd_o}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/streamctl_done_i}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/streamctl_valid_i}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/streamctl_ready_o}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/curr_cfg_d}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/curr_cfg_q}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/curr_cfg_idx_d}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/frep_cnt_d}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/frep_cnt_q}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/inst_last}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/rpt_last}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/seq_done}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/loop_end}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/seq_out_ready}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/seq_out_valid}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/seq_qdata_op}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/seq_qdata_argc}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/core_rb_valid}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/core_rb_ready}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/core_rpt_valid}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/core_direct_valid}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/core_direct_ready}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/rd_pointer_d}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/rd_pointer_q}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/mem_d}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/mem_q}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/rb_wr_pointer}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/rb_rd_data}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/rb_wr_data}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/rb_wr_en}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/rb_full}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/rb_empty}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/rb_contains_instructions}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/rb_contains_no_instructions}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/wr_pointer_d}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/wr_pointer_q}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/rpt_cnt_d}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/inst_cnt_d}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/stagger_cnt_q}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/stagger_cnt_d}
add wave -noupdate -color Orange {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/oup_qdata_op_o}
add wave -noupdate -color Orange {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/oup_qvalid_o}
add wave -noupdate -color Orange {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/oup_qready_i}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/curr_cfg_idx_q}
add wave -noupdate -expand {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/rpt_cnt_q}
add wave -noupdate -expand {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/inst_cnt_q}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/seq_next}
TreeUpdate [SetDefaultTree]
WaveRestoreCursors {{Cursor 1} {658865 ps} 0}
quietly wave cursor active 1
configure wave -namecolwidth 150
configure wave -valuecolwidth 100
configure wave -justifyvalue left
configure wave -signalnamewidth 1
configure wave -snapdistance 10
configure wave -datasetprefix 0
configure wave -rowmargin 4
configure wave -childrowmargin 2
configure wave -gridoffset 0
configure wave -gridperiod 1
configure wave -griddelta 40
configure wave -timeline 0
configure wave -timelineunits ns
update
WaveRestoreZoom {649569 ps} {668299 ps}
