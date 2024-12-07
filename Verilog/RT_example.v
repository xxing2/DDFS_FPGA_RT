////////////////////////////////////////////////////////////////////////////////
// This Verilog code is an example implementation of the Recursive Trigonometry
// (RT) technique applied to Direct Digital Frequency Synthesizers (DDFS) on FPGA.
// It corresponds to the approach presented in the following publication:
//
// A Recursive Trigonometric Technique for Direct Digital Frequency Synthesizer Implementation
// by Xing Xing, William Melek, and Wilson Wang
// Department of Mechanical and Mechatronics Engineering, University of Waterloo, Waterloo, ON N2L 3G1, Canada
// Department of Mechanical and Mechatronics Engineering, Lakehead University, Thunder Bay, ON P7B 5E1, Canada
// Electronics 2024, 13(23), 4762; https://doi.org/10.3390/electronics13234762
// Submission received: 22 October 2024 / Revised: 25 November 2024
// Accepted: 1 December 2024 / Published: 2 December 2024
//
// Original code authors: Xing Xing, William Melek, and Wilson Wang
// Code optimization for FPGA implementation: Yanan Gao (NVIDIA Senior Engineer)
//
// Usage or reference of this code should cite the above paper.
////////////////////////////////////////////////////////////////////////////////


module cos #(
    parameter   DATA_W      = 8
)
(
    input   wire                        rst,
    input   wire                        clk,
    output  wire                        dac_mode,
    output  wire                        dac_clka,
    output  reg     [DATA_W - 1 : 0]    dac_da,
    output  wire                        dac_wra,
    output  wire                        dac_sleep
);

reg signed        [31:0]   cos;
reg signed        [7:0]   cos1;
reg signed        [7:0]   cos2;
reg signed        [7:0]   cos3;


    //
    // counter
    //
    localparam  COUNT_MAX   = 72;
    localparam  COUNT_W     = $clog2(COUNT_MAX);

    reg     [COUNT_W - 1 : 0]   r_counter;
    wire    [COUNT_W - 1 : 0]   w_counter_next;
    wire                        w_counter_tick;

    always @ (posedge clk, posedge rst)
        if (~rst)
            r_counter <= {COUNT_W{1'b0}};
        else
            r_counter <= w_counter_next;

    assign w_counter_tick = (r_counter == (COUNT_MAX[COUNT_W - 1 : 0] - 1'b1));
    assign w_counter_next = w_counter_tick ? {COUNT_W{1'b0}} : r_counter + 1'b1;


    //
    //
    //
    localparam  REG_W   = 32;

    localparam  COS_0   = 32'h40000000;
    localparam  COS_A   = 32'h3fc1a768;
    localparam  C       = 32'h7f834ed0;

    reg signed  [REG_W - 1 : 0]     r_first;
    reg signed  [REG_W - 1 : 0]     r_second;
    reg signed  [REG_W - 1 : 0]     w_result;
    reg signed  [2 * REG_W - 1 : 0] w_product;
    reg                             w_load;

    always @ (posedge clk, posedge rst)
        if (~rst) begin
            r_first  <= COS_0;
            r_second <= COS_A;
        end else begin
            if (w_load) begin
                r_first  <= COS_0;
                r_second <= COS_A;
            end else begin
                r_first  <= r_second;
                r_second <= w_result;
            end
        end

    always @(*) begin
        w_product = r_second * $signed(C);

        w_result = w_product[2 * REG_W - 1 : 30] - r_first;

        w_load = w_counter_tick;

        cos = r_first;
    end

//////////////////////////signal dac_da gets from 255-data/////////////////////////
always @ (posedge clk or negedge rst)
begin
    if(!rst)
        begin
            cos1   <= 1'b0;
        end
    else
        begin
            cos1   <= cos[31:24];
        end
end

always @ (posedge clk or negedge rst)
begin
    if(!rst)
    begin
        cos2 <= 1'b0;
    end
    else
    begin 
        cos2 <= (cos1>>7);
    end
end

always @ (*)
begin
    if(cos2 == 7'b0)
        cos3 = cos1 + 7'd255;
    else
        cos3 = cos1 - 7'd255;
end
///////////////////////////////////////////////////////////////////////////////////
always  @(posedge clk or negedge rst)begin
    if(rst==1'b0)begin
        dac_da <= 0;
    end
    else begin
        dac_da <= 255-cos3;//255-cos3;
    end
end

assign dac_sleep = 0      ;
assign dac_wra = dac_clka ;
assign dac_clka = ~clk    ;
assign dac_mode = 1       ;


endmodule

