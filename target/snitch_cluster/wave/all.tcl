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
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/oup_qdata_op_o}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/oup_qdata_arga_o}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/oup_qdata_argb_o}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/oup_qdata_argc_o}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/oup_qdata_repd_o}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/oup_qvalid_o}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/oup_qready_i}
add wave -noupdate -color {Cornflower Blue} -expand -subitemconfig {{/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/frep_cfg_q[3]} {-color {Cornflower Blue} -height 16} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/frep_cfg_q[2]} {-color {Cornflower Blue} -height 16} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/frep_cfg_q[1]} {-color {Cornflower Blue} -height 16} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/frep_cfg_q[0]} {-color {Cornflower Blue} -height 16}} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/frep_cfg_q}
add wave -noupdate -color Orange -radix binary {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/last_inst}
add wave -noupdate -color Orange -radix binary {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/last_iter}
add wave -noupdate -color {Cornflower Blue} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/frep_idx_q}
add wave -noupdate -color {Cornflower Blue} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/frep_cnt_q}
add wave -noupdate -radix binary {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/inst_starts_loop}
add wave -noupdate -radix binary {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/inst_starts_mask}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/starting_frep_idx}
add wave -noupdate -color Orange {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/seq_next}
add wave -noupdate -color Orange {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/frep_active_d}
add wave -noupdate -color Orange {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/frep_active_q}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/seq_done}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/seq_out_ready}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/seq_out_valid}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/seq_qdata_op}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/seq_qdata_argc}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/frep_is_nested}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/core_rb_valid}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/core_rb_ready}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/core_rpt_valid}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/core_direct_valid}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/core_direct_ready}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/rd_pointer_d}
add wave -noupdate -color {Cornflower Blue} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/rd_pointer_q}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/rb_empty}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/stagger_cnt_q}
add wave -noupdate {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/stagger_cnt_d}
add wave -noupdate -expand -group {ring buffer} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/i_ring_buffer/clk_i}
add wave -noupdate -expand -group {ring buffer} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/i_ring_buffer/rst_ni}
add wave -noupdate -expand -group {ring buffer} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/i_ring_buffer/wvalid_i}
add wave -noupdate -expand -group {ring buffer} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/i_ring_buffer/wready_o}
add wave -noupdate -expand -group {ring buffer} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/i_ring_buffer/wdata_i}
add wave -noupdate -expand -group {ring buffer} -color Orange {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/i_ring_buffer/raddr_i}
add wave -noupdate -expand -group {ring buffer} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/i_ring_buffer/rdata_o}
add wave -noupdate -expand -group {ring buffer} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/i_ring_buffer/advance_i}
add wave -noupdate -expand -group {ring buffer} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/i_ring_buffer/step_i}
add wave -noupdate -expand -group {ring buffer} -color Orange {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/i_ring_buffer/rptr_o}
add wave -noupdate -expand -group {ring buffer} -color Orange {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/i_ring_buffer/wptr_o}
add wave -noupdate -expand -group {ring buffer} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/i_ring_buffer/full_o}
add wave -noupdate -expand -group {ring buffer} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/i_ring_buffer/empty_o}
add wave -noupdate -expand -group {ring buffer} -color {Cornflower Blue} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/i_ring_buffer/rptr_q}
add wave -noupdate -expand -group {ring buffer} -color {Cornflower Blue} {/tb_bin/i_dut/i_snitch_cluster/i_cluster/gen_core[0]/i_snitch_cc/gen_fpu/i_snitch_fp_ss/gen_fpu_sequencer/i_snitch_fpu_sequencer/i_ring_buffer/wptr_q}
TreeUpdate [SetDefaultTree]
WaveRestoreCursors {{FMUL 0} {891000 ps} 1} {FREP {903000 ps} 1} {{Out of bounds read} {912000 ps} 1}
quietly wave cursor active 3
configure wave -namecolwidth 210
configure wave -valuecolwidth 233
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
WaveRestoreZoom {892191 ps} {925751 ps}
