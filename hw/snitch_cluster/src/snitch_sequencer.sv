// Copyright 2020 ETH Zurich and University of Bologna.
// Solderpad Hardware License, Version 0.51, see LICENSE for details.
// SPDX-License-Identifier: SHL-0.51

`include "common_cells/registers.svh"
`include "common_cells/assertions.svh"

/// Description: Filters FPU repetition instructions
module snitch_sequencer import snitch_pkg::*; #(
    parameter int unsigned AddrWidth = 0,
    parameter int unsigned DataWidth = 0,
    parameter int unsigned Depth = 32,
    parameter acc_addr_e DstAddr = FP_SS,
    /// Derived parameter *Do not override*
    localparam type addr_t = logic [AddrWidth-1:0],
    localparam type data_t = logic [DataWidth-1:0]
) (
    input  logic                             clk_i,
    input  logic                             rst_i,
    // pragma translate_off
    output fpu_sequencer_trace_port_t        trace_port_o,
    // pragma translate_on
    input  acc_addr_e                        inp_qaddr_i,
    input  logic                      [ 4:0] inp_qid_i,
    input  logic                      [31:0] inp_qdata_op_i,  // RISC-V instruction
    input  data_t                            inp_qdata_arga_i,
    input  data_t                            inp_qdata_argb_i,
    input  addr_t                            inp_qdata_argc_i,
    input  logic                             inp_qvalid_i,
    output logic                             inp_qready_o,
    output acc_addr_e                        oup_qaddr_o,
    output logic                      [ 4:0] oup_qid_o,
    output logic                      [31:0] oup_qdata_op_o,  // RISC-V instruction
    output data_t                            oup_qdata_arga_o,
    output data_t                            oup_qdata_argb_o,
    output addr_t                            oup_qdata_argc_o,
    output logic                             oup_qdata_repd_o,  // Whether this is a repeated issue
    output logic                             oup_qvalid_o,
    input  logic                             oup_qready_i,
    // SSR stream control interface
    input  logic                             streamctl_done_i,
    input  logic                             streamctl_valid_i,
    output logic                             streamctl_ready_o
);

  // TODO: uniformize name variants {frep, seq}, {rpt, iter}, {direct, bypass}
  // TODO: check all always_comb blocks have meaningful names

  localparam int RptBits = 16;
  localparam int FrepDim = 4;

  localparam int unsigned DepthBits = $clog2(Depth);
  // TODO: check which types can be loosened to use LoopIdxBits
  localparam int unsigned LoopCntBits = $clog2(FrepDim + 1);
  localparam int unsigned LoopIdxBits = $clog2(FrepDim);

  /////////////
  // Decoder //
  /////////////

  typedef enum logic [1:0] {
    Frep = 2,
    Bypass = 1,
    Buffer = 0
  } seq_path_e;

  seq_path_e inst_path_select;

  always_comb begin
    inst_path_select = Buffer;

    unique casez (inp_qdata_op_i)
      riscv_instr::FREP_O,
      riscv_instr::FREP_I, riscv_instr::IREP: begin
        inst_path_select = Frep;
      end

      // Instructions which explicitly sync between int and float
      // pipeline are not supported within an FREP and need to
      // bypass the ring buffer.

      // float to int
      riscv_instr::FLE_S,
      riscv_instr::FLT_S,
      riscv_instr::FEQ_S,
      riscv_instr::FCLASS_S,
      riscv_instr::FCVT_W_S,
      riscv_instr::FCVT_WU_S,
      riscv_instr::FMV_X_W,
      riscv_instr::VFEQ_S,
      riscv_instr::VFEQ_R_S,
      riscv_instr::VFNE_S,
      riscv_instr::VFNE_R_S,
      riscv_instr::VFLT_S,
      riscv_instr::VFLT_R_S,
      riscv_instr::VFGE_S,
      riscv_instr::VFGE_R_S,
      riscv_instr::VFLE_S,
      riscv_instr::VFLE_R_S,
      riscv_instr::VFGT_S,
      riscv_instr::VFGT_R_S,
      riscv_instr::VFCLASS_S,
      riscv_instr::FLE_D,
      riscv_instr::FLT_D,
      riscv_instr::FEQ_D,
      riscv_instr::FCLASS_D,
      riscv_instr::FCVT_W_D,
      riscv_instr::FCVT_WU_D,
      riscv_instr::FMV_X_D,
      riscv_instr::FLE_H,
      riscv_instr::FLT_H,
      riscv_instr::FEQ_H,
      riscv_instr::FCLASS_H,
      riscv_instr::FCVT_W_H,
      riscv_instr::FCVT_WU_H,
      riscv_instr::FMV_X_H,
      riscv_instr::VFEQ_H,
      riscv_instr::VFEQ_R_H,
      riscv_instr::VFNE_H,
      riscv_instr::VFNE_R_H,
      riscv_instr::VFLT_H,
      riscv_instr::VFLT_R_H,
      riscv_instr::VFGE_H,
      riscv_instr::VFGE_R_H,
      riscv_instr::VFLE_H,
      riscv_instr::VFLE_R_H,
      riscv_instr::VFGT_H,
      riscv_instr::VFGT_R_H,
      riscv_instr::VFCLASS_H,
      riscv_instr::VFMV_X_H,
      riscv_instr::VFCVT_X_H,
      riscv_instr::VFCVT_XU_H,
      riscv_instr::FLE_AH,
      riscv_instr::FLT_AH,
      riscv_instr::FEQ_AH,
      riscv_instr::FCLASS_AH,
      riscv_instr::FMV_X_AH,
      riscv_instr::VFEQ_AH,
      riscv_instr::VFEQ_R_AH,
      riscv_instr::VFNE_AH,
      riscv_instr::VFNE_R_AH,
      riscv_instr::VFLT_AH,
      riscv_instr::VFLT_R_AH,
      riscv_instr::VFGE_AH,
      riscv_instr::VFGE_R_AH,
      riscv_instr::VFLE_AH,
      riscv_instr::VFLE_R_AH,
      riscv_instr::VFGT_AH,
      riscv_instr::VFGT_R_AH,
      riscv_instr::VFCLASS_AH,
      riscv_instr::VFMV_X_AH,
      riscv_instr::VFCVT_X_AH,
      riscv_instr::VFCVT_XU_AH,
      riscv_instr::FLE_B,
      riscv_instr::FLT_B,
      riscv_instr::FEQ_B,
      riscv_instr::FCLASS_B,
      riscv_instr::FCVT_W_B,
      riscv_instr::FCVT_WU_B,
      riscv_instr::FMV_X_B,
      riscv_instr::VFEQ_B,
      riscv_instr::VFEQ_R_B,
      riscv_instr::VFNE_B,
      riscv_instr::VFNE_R_B,
      riscv_instr::VFLT_B,
      riscv_instr::VFLT_R_B,
      riscv_instr::VFGE_B,
      riscv_instr::VFGE_R_B,
      riscv_instr::VFLE_B,
      riscv_instr::VFLE_R_B,
      riscv_instr::VFGT_B,
      riscv_instr::VFGT_R_B,
      riscv_instr::VFCLASS_B,
      riscv_instr::VFMV_X_B,
      riscv_instr::VFCVT_X_B,
      riscv_instr::VFCVT_XU_B,

      // int to float
      riscv_instr::FMV_W_X,
      riscv_instr::FCVT_S_W,
      riscv_instr::FCVT_S_WU,
      riscv_instr::FCVT_D_W,
      riscv_instr::FCVT_D_WU,
      riscv_instr::FMV_H_X,
      riscv_instr::FCVT_H_W,
      riscv_instr::FCVT_H_WU,
      riscv_instr::VFMV_H_X,
      riscv_instr::VFCVT_H_X,
      riscv_instr::VFCVT_H_XU,
      riscv_instr::FMV_AH_X,
      riscv_instr::VFMV_AH_X,
      riscv_instr::VFCVT_AH_X,
      riscv_instr::VFCVT_AH_XU,
      riscv_instr::FMV_B_X,
      riscv_instr::FCVT_B_W,
      riscv_instr::FCVT_B_WU,
      riscv_instr::VFMV_B_X,
      riscv_instr::VFCVT_B_X,
      riscv_instr::VFCVT_B_XU,
      // riscv_instr::IMV_X_W,
      // riscv_instr::IMV_W_X,

      // CSR accesses
      riscv_instr::CSRRW,
      riscv_instr::CSRRS,
      riscv_instr::CSRRC,
      riscv_instr::CSRRWI,
      riscv_instr::CSRRSI,
      riscv_instr::CSRRCI: begin
        inst_path_select = Bypass;
      end

      // All other instructions have to be buffered
      default: begin
        inst_path_select = Buffer;
      end
    endcase
  end

  /////////////////
  // Input demux //
  /////////////////

  logic core_rb_valid, core_rb_ready;
  logic core_rpt_valid, core_rpt_ready;
  logic core_direct_valid, core_direct_ready;

  stream_demux #(
    .N_OUP(3)
  ) i_input_demux (
    .inp_valid_i(inp_qvalid_i),
    .inp_ready_o(inp_qready_o),
    .oup_sel_i(inst_path_select),
    .oup_valid_o({core_rpt_valid, core_direct_valid, core_rb_valid}),
    .oup_ready_i({core_rpt_ready, core_direct_ready, core_rb_ready})
  );

  ////////////////
  // Loop state //
  ////////////////

  logic [FrepDim-1:0] last_iter;
  logic [FrepDim-1:0] last_inst;

  logic [FrepDim-1:0][RptBits-1:0] rpt_cnt;

  // TODO
  logic [2:0] stagger_cnt_q, stagger_cnt_d;
  `FFAR(stagger_cnt_q, stagger_cnt_d, '0, clk_i, rst_i)

  /////////////////////
  // Loop nest state //
  /////////////////////

  typedef struct packed {
    logic is_streamctl;
    logic is_outer;
    logic [DepthBits-1:0] max_inst;
    logic [RptBits-1:0] max_rpt;
    logic [2:0] stagger_max;
    logic [3:0] stagger_mask;  // one-hot stagger mask
    logic [DepthBits:0] base_pointer;  // loop base pointer where the config starts
  } seq_cfg_t;

  seq_cfg_t [FrepDim-1:0] frep_cfg_d, frep_cfg_q;
  logic frep_active_d, frep_active_q;

  logic [LoopCntBits-1:0] frep_idx_d, frep_idx_q;
  logic [LoopCntBits-1:0] frep_cnt_d, frep_cnt_q;

  `FFAR(frep_cfg_q, frep_cfg_d, '0, clk_i, rst_i);
  `FFAR(frep_active_q, frep_active_d, '0, clk_i, rst_i);
  `FFAR(frep_idx_q, frep_idx_d, '0, clk_i, rst_i);
  `FFAR(frep_cnt_q, frep_cnt_d, '0, clk_i, rst_i);

  logic [FrepDim-1:0][DepthBits-1:0] loop_end_pointer;
  // TODO Do we need both loop_end and seq_done?
  logic seq_next, seq_done, loop_end;

  logic seq_out_ready, seq_out_valid;
  logic [31:0] seq_qdata_op;
  logic [AddrWidth-1:0] seq_qdata_argc;

  /////////////////
  // Ring buffer //
  /////////////////

  typedef struct packed {
    logic [31:0] qdata_op;
    addr_t qdata_argc;
  } seq_entry_t;

  seq_entry_t rb_rdata;
  seq_entry_t rb_wdata;
  logic rb_empty;
  logic [DepthBits-1:0] rb_wptr;
  logic rb_advance;
  // TODO comment width
  logic [$clog2(Depth+1)-1:0] rb_step;

  // TODO does this belong here
  logic [DepthBits-1:0] rd_pointer_d, rd_pointer_q;
  `FFAR(rd_pointer_q, rd_pointer_d, '0, clk_i, rst_i)

  assign rb_wdata = '{
    qdata_op: inp_qdata_op_i,
    qdata_argc: inp_qdata_argc_i
  };

  assign rb_advance = frep_active_q ? loop_end : seq_next;
  assign rb_step = frep_active_q ? frep_cfg_q[0].max_inst + 1 : 1;

  ring_buffer #(
    .Depth(Depth),
    .data_t(seq_entry_t)
  ) i_ring_buffer (
    .clk_i(clk_i),
    .rst_ni(~rst_i),
    .wvalid_i(core_rb_valid),
    .wready_o(core_rb_ready),
    .wdata_i(rb_wdata),
    .raddr_i(rd_pointer_q),
    .rdata_o(rb_rdata),
    .advance_i(rb_advance),
    .step_i(rb_step),
    .wptr_o(rb_wptr),
    .rptr_o(),
    .full_o(),
    .empty_o(rb_empty)
  );

  // seq_entry_t [Depth-1:0] mem_d, mem_q;
  // logic [DepthBits:0] wr_pointer_d, wr_pointer_q;

  // always_comb begin : ringbuffer
  //   mem_d = mem_q;
  //   if (rb_wr_en) mem_d[wr_pointer_q[DepthBits-1:0]] = rb_wdata;
  //   rb_rdata = mem_q[rd_pointer_q[DepthBits-1:0]];
  //   rb_full = (rd_pointer_q ^ wr_pointer_q) == 1 << DepthBits;
  //   rb_empty = (rd_pointer_q ^ wr_pointer_q) == '0;
  // end

  // `FFAR(mem_q, mem_d, '0, clk_i, rst_i)
  // `FFAR(wr_pointer_q, wr_pointer_d, '0, clk_i, rst_i)

  // always_comb begin
  //   wr_pointer_d = wr_pointer_q;

  //   rb_wr_en = 0;
  //   core_rb_ready = !rb_full;

  //   // write logic
  //   if (core_rb_valid && core_rb_ready) begin
  //     wr_pointer_d = wr_pointer_q + 1;
  //     rb_wr_en = 1;
  //   end
  // end

  /////////////////////////////////
  // Loop state transition logic //
  /////////////////////////////////

  logic [FrepDim-1:0] last_iter_inner_loops;

  for (genvar i = 0; i < FrepDim; i++) begin : gen_fsm_loop

    logic incr_inst;
    logic incr_iter;

    logic [FrepDim-1:0] loop_mask;

    // // Mask over present and outer loops, and inner inactive loops
    // // (i.e. loops nested within the currently active loop).
    // mask #(
    //   .Width(FrepDim),
    //   .Mode(1)
    // ) i_mask (
    //   .lsb_i(i + 1),
    //   .msb_i(frep_idx_q),
    //   .mask_o(loop_mask)
    // );

    // Instructions in inner loop bodies may be repeated multiple times, but must
    // be counted only once by the outer loop instruction counter. Specifically,
    // we count them when they are last issued, i.e. when all loops between the
    // one containing the current instruction (frep_idx_q) and the outer loop (i)
    // are in their last iteration.
    logic [FrepDim-1:0] outer_loops_mask, inner_loops_mask;
    assign outer_loops_mask = (1 << (i + 1)) - 1;  // Ignore present and outer loops
    assign inner_loops_mask = ~((1 << (frep_idx_q + 1)) - 1);  // Ignore inner, inactive loops
    // assign last_iter_inner_loops[i] = &(last_iter | loop_mask);
    assign last_iter_inner_loops[i] = &(last_iter | outer_loops_mask | inner_loops_mask);
    assign incr_inst = frep_active_q && seq_next
      && (((i == frep_idx_q) && (i < frep_cnt_q))
      || ((frep_idx_q > i) && last_iter_inner_loops[i]));

    trip_counter #(
      .WIDTH(DepthBits)
    ) i_inst_counter (
      .clk_i(clk_i),
      .rst_ni(~rst_i),
      .en_i(incr_inst),
      .delta_i(DepthBits'(1)),
      .bound_i(frep_cfg_q[i].max_inst),
      .q_o(),
      .last_o(last_inst[i]),
      .trip_o(incr_iter)
    );

    trip_counter #(
      .WIDTH(RptBits)
    ) i_rpt_counter (
      .clk_i(clk_i),
      .rst_ni(~rst_i),
      .en_i(incr_iter),
      .delta_i(RptBits'(1)),
      .bound_i(frep_cfg_q[i].max_rpt),
      .q_o(rpt_cnt[i]),
      .last_o(last_iter[i]),
      .trip_o()
    );

    // TODO does this belong here?
    assign loop_end_pointer[i] = frep_cfg_q[i].base_pointer + frep_cfg_q[i].max_inst;

  end

  //////////////////////////////////////
  // Loop nest state transition logic //
  //////////////////////////////////////

  logic frep_is_nested;

  // FREP configuration and counter update
  always_comb begin
    frep_cfg_d = frep_cfg_q;
    frep_cnt_d = frep_cnt_q;

    // TODO does this belong here?
    // FREPs can only be accepted if no loop nest is currently configured, or if it
    // is nested within the current loop nest.
    if (frep_cfg_q[0].base_pointer <= loop_end_pointer[0]) begin
      frep_is_nested = (rb_wptr >= frep_cfg_q[0].base_pointer) && (rb_wptr <= loop_end_pointer[0]);
    end else begin
      frep_is_nested = (rb_wptr >= frep_cfg_q[0].base_pointer) || (rb_wptr <= loop_end_pointer[0]);
    end
    core_rpt_ready = (frep_cnt_q == 0) || frep_is_nested;

    if (core_rpt_valid && core_rpt_ready) begin
      frep_cfg_d[frep_cnt_q].is_streamctl = inp_qdata_op_i[31];
      frep_cfg_d[frep_cnt_q].is_outer = inp_qdata_op_i[7];
      frep_cfg_d[frep_cnt_q].max_inst = inp_qdata_op_i[20+:DepthBits];
      frep_cfg_d[frep_cnt_q].stagger_mask = inp_qdata_op_i[11:8];
      frep_cfg_d[frep_cnt_q].stagger_max = inp_qdata_op_i[14:12];
      frep_cfg_d[frep_cnt_q].max_rpt = inp_qdata_arga_i[RptBits-1:0];
      frep_cfg_d[frep_cnt_q].base_pointer = rb_wptr;

      frep_cnt_d = frep_cnt_q + 1;
    end 

    if (loop_end) begin
      frep_cnt_d = 0;
    end
  end

  logic [FrepDim-1:0] inst_starts_loop;

  for (genvar i = 0; i < FrepDim; i++) begin : gen_inst_starts_loop
    assign inst_starts_loop[i] = (rd_pointer_d == frep_cfg_d[i].base_pointer);
  end

  // // Mask over 
  // mask #(
  //   .Width(FrepDim),
  //   .Mode(1)
  // ) i_mask (
  //   .lsb_i(i + 1),      // Ignore present and outer loops
  //   .msb_i(frep_idx_q), // Ignore inner, inactive loops
  //   .mask_o(mask)
  // );

  logic [FrepDim-1:0] inst_starts_lmask, inst_starts_hmask, inst_starts_mask;
  logic no_loop_starts;
  logic [LoopIdxBits-1:0] lzc_cnt;
  logic [LoopIdxBits-1:0] starting_frep_idx;

  // Mask to select only the loops which need to be considered
  // to find the innermost loop starting at the next instruction,
  // That is all loops i, with frep_idx_q < i < frep_cnt_d.
  assign inst_starts_lmask = (1 << (frep_idx_q + 1)) - 1;
  assign inst_starts_hmask = ~((1 << (frep_cnt_d - 1 + 1)) - 1);
  assign inst_starts_mask = ~(inst_starts_lmask | inst_starts_hmask);

  // Find the innermost loop starting at the next instruction.
  lzc #(
    .WIDTH(FrepDim),
    .MODE(1)
  ) i_loop_start_lzc (
    .in_i(inst_starts_mask & inst_starts_loop),
    .cnt_o(lzc_cnt),
    .empty_o(no_loop_starts)
  );
  assign starting_frep_idx = FrepDim - lzc_cnt - 1;

  logic [LoopIdxBits-1:0] non_ending_loops_cnt;
  logic [LoopIdxBits-1:0] outermost_non_ending_loop;
  logic [FrepDim-1:0] loop_ends;
  logic [FrepDim-1:0] loop_active;
  logic no_loop_ends;

  assign loop_ends = last_inst & last_iter & last_iter_inner_loops;
  assign loop_active = (1 << (frep_idx_q + 1)) - 1;

  // Compute the innermost active loop which does not end with the current instruction, using a trailing zero counter.
  lzc #(
    .WIDTH(FrepDim),
    .MODE(0)
  ) i_loop_end_tzc (
    .in_i(loop_ends | ~loop_active),
    .cnt_o(non_ending_loops_cnt),
    .empty_o(no_loop_ends)
  );
  assign outermost_non_ending_loop = no_loop_ends ? non_ending_loops_cnt : non_ending_loops_cnt - 1;

  // FREP index update
  always_comb begin : sequence_logic
    frep_idx_d = frep_idx_q;
    rd_pointer_d = rd_pointer_q;
    frep_active_d = frep_active_q;
    loop_end = 1'b0;

    // Update read pointer into the ring buffer
    if (seq_next) begin
      rd_pointer_d = rd_pointer_q + 1;
    end

    // Update frep active flag
    if (frep_cnt_d > 0 && inst_starts_loop[0]) begin
      frep_active_d = 1'b1;
    end

    // As we complete an inner loop we must update the frep_idx to the outermost
    // non-complete loop
    if (frep_active_q && seq_next && frep_cfg_q[frep_idx_q].is_outer) begin
      // Reset loop nest if all loops are at the end
      if (non_ending_loops_cnt == 0) begin
        frep_idx_d = '0;
        loop_end = 1'b1;
        frep_active_d = 1'b0;
      end else begin
        // Otherwise move to the outermost non-ending loop
        frep_idx_d = outermost_non_ending_loop;
        // If this loop is also at the last instruction, we must reset the
        // read pointer to the base of the loop.
        if (last_inst[outermost_non_ending_loop]) begin
          rd_pointer_d = frep_cfg_q[outermost_non_ending_loop].base_pointer;
        end
      end
    end

    // As we move to the next instruction we must update the frep_idx
    // to the innermost starting loop.
    if (frep_active_q && !no_loop_starts) begin
      frep_idx_d = starting_frep_idx;
    end
  end

  ////////////////////////////
  // Loop nest output logic //
  ////////////////////////////

  assign seq_next = seq_out_valid & seq_out_ready;

  // TODO understand
  always_comb begin : proc_streamctl
    seq_out_valid     = !rb_empty;
    seq_done          = 1'b0;
    streamctl_ready_o = 1'b0;
    if ((frep_cnt_q > 0) && frep_cfg_q[frep_idx_q].is_outer && frep_cfg_q[frep_idx_q].is_streamctl) begin
      seq_out_valid     = !rb_empty && streamctl_valid_i && !streamctl_done_i;
      seq_done          = !rb_empty && streamctl_valid_i && streamctl_done_i;
      streamctl_ready_o = (!rb_empty && seq_out_ready) || seq_done;
    end
  end

  // Compose offloading instruction e.g. staggering.
  always_comb begin
    seq_qdata_op   = rb_rdata.qdata_op;
    seq_qdata_argc = rb_rdata.qdata_argc;
    // TODO
    // if (frep_cfg_q[frep_idx_q].stagger_mask[0]) seq_qdata_op[11:7] += stagger_cnt_q;
    // if (frep_cfg_q[frep_idx_q].stagger_mask[1]) seq_qdata_op[19:15] += stagger_cnt_q;
    // if (frep_cfg_q[frep_idx_q].stagger_mask[2]) seq_qdata_op[24:20] += stagger_cnt_q;
    // if (frep_cfg_q[frep_idx_q].stagger_mask[3]) seq_qdata_op[31:27] += stagger_cnt_q;
  end

  ////////////////
  // Output mux //
  ////////////////

  typedef struct packed {
    acc_addr_e   qaddr;
    logic  [4:0] qid;
    logic [31:0] qdata_op;  // RISC-V instruction
    data_t       qdata_arga;
    data_t       qdata_argb;
    addr_t       qdata_argc;
    logic        qdata_repd;
  } seq_data_t;

  seq_data_t core_direct_data, seq_out_data, oup_data;

  assign core_direct_data = '{
    qaddr:      inp_qaddr_i,
    qid:        inp_qid_i,
    qdata_op:   inp_qdata_op_i,
    qdata_arga: inp_qdata_arga_i,
    qdata_argb: inp_qdata_argb_i,
    qdata_argc: inp_qdata_argc_i,
    qdata_repd: 1'b0
  };

  assign seq_out_data = '{
    qaddr:      DstAddr,
    qid:        '0,
    qdata_op:   seq_qdata_op,
    qdata_arga: '0,
    qdata_argb: '0,
    qdata_argc: $unsigned(seq_qdata_argc),
    // If this repeats a previously issued instruction, communicate this
    // to subsystem (e.g. for single issuing of CAQ responses).
    qdata_repd: (rpt_cnt[frep_idx_q] != 0)
  };

  // Select bypass path iff ring buffer is empty
  stream_mux #(
    .DATA_T(seq_data_t),
    .N_INP(2)
  ) i_output_mux (
    .inp_data_i({core_direct_data, seq_out_data}),
    .inp_valid_i({core_direct_valid, seq_out_valid}),
    .inp_ready_o({core_direct_ready, seq_out_ready}),
    .inp_sel_i(rb_empty),
    .oup_data_o(oup_data),
    .oup_valid_o(oup_qvalid_o),
    .oup_ready_i(oup_qready_i)
  );

  assign oup_qaddr_o      = oup_data.qaddr;
  assign oup_qid_o        = oup_data.qid;
  assign oup_qdata_op_o   = oup_data.qdata_op;
  assign oup_qdata_arga_o = oup_data.qdata_arga;
  assign oup_qdata_argb_o = oup_data.qdata_argb;
  assign oup_qdata_argc_o = oup_data.qdata_argc;
  assign oup_qdata_repd_o = oup_data.qdata_repd;

  ////////////
  // Tracer //
  ////////////

  // pragma translate_off
  assign trace_port_o.source    = snitch_pkg::SrcFpuSeq;
  assign trace_port_o.cbuf_push = core_rpt_valid && core_rpt_ready;
  assign trace_port_o.is_outer  = frep_cfg_d[frep_cnt_q].is_outer;
  assign trace_port_o.max_inst  = frep_cfg_d[frep_cnt_q].max_inst;
  assign trace_port_o.max_rpt   = frep_cfg_d[frep_cnt_q].max_rpt;
  assign trace_port_o.stg_max   = frep_cfg_d[frep_cnt_q].stagger_max;
  assign trace_port_o.stg_mask  = frep_cfg_d[frep_cnt_q].stagger_mask;
  // pragma translate_on

  ////////////////
  // Assertions //
  ////////////////

  // TODO add assertion to check that frep_cnt never overflows, i.e. we
  // don't have more nested FREPs than supported

  // TODO add assertion to check that, if a loop is active, and not at the last instruction
  // the outer loops must also not be at the last instruction

  // Ensure that `max_inst` bits fit into assigned slot
  `ASSERT_INIT(CheckMaxInstFieldWidth, DepthBits < 11);

endmodule
