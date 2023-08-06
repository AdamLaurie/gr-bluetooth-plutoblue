PLUTO_IP="$1"
scp system_top.bit.bin root@"$PLUTO_IP":/lib/firmware/
scp loadplutobluez_runme_on_pluto.txt root@"$PLUTO_IP":/configurepluto.sh
ssh root@"$PLUTO_IP" "ash /configurepluto.sh"
IIOD_REMOTE=ip:"$1" iio_readdev -T 100000 -b 4000000 -s 128000000 cf-ad9361-lpc voltage0 voltage1 | ../plutoblueconverter/target/release/plutoblueconverter  | ./apps/btrx -s -i - -f 2440000000.0 -S -r 84000000 -t 10
