function socket_sender_udp(out)
u = udp('255.255.255.255', 8053, 'LocalPort', 3112);
fopen(u);
data = int2str(out);
fwrite(u, data, 'char');
fclose(u);
end