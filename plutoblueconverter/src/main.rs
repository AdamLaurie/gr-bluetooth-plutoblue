use arrayfire::*;
use num_complex::*;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use std::io::BufWriter;
const buf_size: u64= 4*1024*1024;
const buf_size_m1: u64=buf_size-1;
fn main() {
    let mut reader = BufReader::with_capacity(128*1024*1024,std::io::stdin());
    let mut writer = BufWriter::with_capacity(128*1024*1024,std::io::stdout());
    let mut buffer:Vec<u8>=vec![0u8;4*1024*1024];
    let mut buffer_count=0;
    let mut sample_arrs=vec![constant!(1u16;2*1024*1024),constant!(2u16;2*1024*1024),constant!(4u16;2*1024*1024),constant!(8u16;2*1024*1024),constant!(16u16;2*1024*1024),constant!(32u16;2*1024*1024),constant!(64u16;2*1024*1024),constant!(128u16;2*1024*1024),constant!(256u16;2*1024*1024),constant!(512u16;2*1024*1024),constant!(1024u16;2*1024*1024),constant!(2048u16;2*1024*1024),constant!(4096u16;2*1024*1024),constant!(8192u16;2*1024*1024),constant!(16384u16;2*1024*1024),constant!(32768u16;2*1024*1024)];
    let mut high_arrs=vec![];
    for i in 0..16 {
        high_arrs.push(constant!(32767i16;2*1024*1024));
    }
    let mut low_arrs=vec![];
    for i in 0..16 {
        low_arrs.push(constant!(-32767i16;2*1024*1024));
    }
    let low_arr=arrayfire::join_many(1i32,low_arrs.iter().map(|a|a).collect());
    let high_arr=arrayfire::join_many(1i32,high_arrs.iter().map(|a|a).collect());
    loop {
        if (buffer_count == 4*1024*1024) {
            let mut gpu_arr=Array::new(&buffer[..],dim4!(4*1024*1024));
            let mut u16_arr=gpu_arr.cast::<u16>();
            let mut u16_upper_arr=gpu_arr.copy();
            let upper = Seq::new(1f64,buf_size_m1 as f64,2f64); 
            let lower = Seq::new(0f64,buf_size_m1 as f64,2f64);

            let mut combined_arr=arrayfire::bitor(&arrayfire::mul(&index(&u16_upper_arr,&[upper]),&constant!(256u16;2*1024*1024),true),&index(&u16_arr,&[lower]),false);
            let mut sample_bit_arr=arrayfire::join_many(1i32,sample_arrs.iter().map(|a| a).collect());
            let mut sample_high=arrayfire::gt(&arrayfire::bitand(&sample_bit_arr,&combined_arr,true),&constant!(0u16;2*1024*1024),true);
            let mut samples=high_arr.copy();           
            arrayfire::replace(&mut samples,&sample_high,&low_arr);
            buffer_count=0;    
            let mut samples_host=vec![0i16;2*1024*1024*16];
            arrayfire::flat(&samples).host(&mut samples_host);
            let mut packed_idx=0;
            while (packed_idx < 1*1024*1024) {
                for samp_idx in 0..16 {
                   writer.write_all(&samples_host[packed_idx*32 + samp_idx].to_le_bytes());
                   writer.write_all(&samples_host[packed_idx*32 + 16 + samp_idx].to_le_bytes());
                }
                packed_idx+=1;
            }
        }
        buffer_count+=reader.read(&mut buffer[buffer_count..]).unwrap();

    }
}
