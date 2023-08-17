`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 08/06/2023 06:18:23 AM
// Design Name: 
// Module Name: sample_quantizer_packer
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module sample_quantizer_packer(
    input aclk,
    input tvalid_s,
    output reg tvalid_m,
    input [15:0] r,
    input [15:0] i,
    output reg [15:0] packed_real,
    output reg [15:0] packed_imag
);
reg [7:0]pack_counter=8'd0;
always @ (posedge aclk) begin
    if (tvalid_s == 1) begin
        packed_real[pack_counter] <= r[15] ? 0 : 1;
        packed_imag[pack_counter] <= i[15] ? 0 : 1;
        if (pack_counter == 15) begin
            pack_counter <= 0;
            tvalid_m <= 1;
        end else begin
            tvalid_m <= 0;
            pack_counter <= pack_counter+1;
        end
    end else begin
        tvalid_m <= 0;
    end
end
endmodule
