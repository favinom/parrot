function faimesh(be,fe)

nameInput='second.txt';

P=load(nameInput);

P=moveNodesPerfect(P,be,fe);

save(nameInput,'P','-ascii');

quit

return



