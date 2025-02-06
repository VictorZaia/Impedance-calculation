# %%
from Config.config import load_config
from Solver.Processor import *
from Solver.Post_Processor import Post_Processor as post
from Solver.Optimization import *
from Solver.Runner import Run
# %% [markdown]
"""
# **Computing and ploting the absorption coefficient for varying parameters**
The impedance calculation will be given by: 

$$ r(w,|v|) = A(w) + B \frac{|v|}{\sigma} + 0.3 * \frac{(1 - \sigma^2)}{\sigma} * M $$

the resistance can be given by a linear part, called A(w) and a non linear part which depends on the acoustic velocity |v|, multiplied by a coefficient B.

$$ r = A(w) + B * \frac{|p|}{\rho  c  \sigma \sqrt{r^2 + (\chi - j \cot{kL})^2}} + 0.3 * \frac{(1 - \sigma^2)}{\sigma} * M $$
"""

# %%
def run():
    print("=========================Importing data=========================")
    runner = Run()
    print("=========================Data imported=========================")

    print("\n=========================Starting impedance study=========================")
    runner.run_impedance_study()
    print("=========================Impedance study finished=========================")

    print("\n=========================Starting optimization study=========================")
    runner.run_optimization()
    print("\n=========================Optimization study finished=========================")

if __name__ == "__main__":
    run()