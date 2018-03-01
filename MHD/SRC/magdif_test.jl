using Plots
pyplot()
#%%
t = ccall((:__magdif_MOD_init, "./libmagdif.so"), Int32, ())
println(t)
#%%
n = unsafe_load(cglobal((:__magdif_MOD_n, "./libmagdif.so"), Int32))
psimin = unsafe_load(cglobal((:__magdif_MOD_psimin, "./libmagdif.so"), Float64))
psimax = unsafe_load(cglobal((:__magdif_MOD_psimax, "./libmagdif.so"), Float64))
nr_max = 74 # TODO: read from config
#%%

psimaxptr = cglobal((:__magdif_MOD_psimax, "./libmagdif.so"), Float64)
q = unsafe_wrap(Array, unsafe_load(cglobal((:__magdif_MOD_q, "./libmagdif.so"),Ptr{Float64})), nr_max)

plot(q)
