# testing exec()

#%% fails
def run_code():
    exec("x = 5")
    print(x)  

run_code()

#%% works

def run_code():
    local_vars = {}
    exec("x = 5", globals(), local_vars)
    print(local_vars["x"])
# %%
run_code()
