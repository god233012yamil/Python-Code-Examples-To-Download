from PySpice.Spice.Netlist import Circuit, SubCircuitFactory
from PySpice.Probe.WaveForm import TransientAnalysis
import numpy as np
import matplotlib.pyplot as plt


class CustomOpAmp(SubCircuitFactory):
    """Custom Op-amp model as a sub circuit"""
    NAME = 'my_op-amp'
    NODES = ('in_plus', 'in_minus', 'pos', 'neg', 'out')

    def __init__(self):
        super().__init__()

        # Define the op-amp parameters.
        A_ol = 1e5       # Open-loop gain (100,000)
        f_unity = 1e6    # Unity-gain bandwidth (1MHz)
        Rf = 1e3         # frequency pole resistor

        # Circuit to simulate the positive supply.
        self.R('pos', 'pos', self.gnd, 1e6)  # 1MΩ

        # Circuit to simulate the negative supply.
        self.R('neg', 'neg', self.gnd, 1e6)  # 1MΩ

        # Define the input impedance.
        self.R('in', 'in_plus', 'in_minus', 1e6)  # 1MΩ

        # Define the op-amp using a BVS (Behavioral Voltage Source).
        # V={A0*V(in+,in-)}
        self.B('B1', 'n1', self.gnd, v='{0}*V(in_plus,in_minus)'.format(A_ol))

        # Calculate the compensation compensation capacitor value to 
        # model finite bandwidth.
        C_comp = 1 / (2 * np.pi * Rf * (f_unity / A_ol))

        # Connect a 1KΩ resistor between node n1 and node n2. 
        self.R('f', 'n1', 'n2', 1e3)  # 1KΩ

        # Connect a C_comp capacitor between node n2 and gnd.
        self.C('f', 'n2', self.gnd, C_comp)

        # Define the output stage using a BVS (Behavioral Voltage Source).
        self.BehavioralSource('B2', 'n3', self.gnd, voltage_expression='V(n2)') 

        # Define the output impedance.
        # Connect a 75Ω resistor between node n3 and the node out. 
        self.R('out', 'n3', 'out', 75)  # 75Ω


def create_circuit() -> Circuit:
    # Create an Inverting op-amp Amplifier circuit.
    circuit = Circuit('Non-inverting Amplifier')

    # Adds the defined CustomOpAmp sub-circuit to the circuit.
    circuit.subcircuit(CustomOpAmp())

    # Connect the my_op-amp pins to circuit nodes.
    circuit.X('1', 'my_op-amp', 
              'in_plus',    # non-inverting input
              'in_minus',   # inverting input
              'pos',        # positive supply
              'neg',        # negative supply
              'out')        # output

    # Connect feedback resistor between the output node (out) and 
    # inverting input node (in_minus).
    circuit.R('1', 'in_minus', 'out', 10e3)
    
    # Connect the input resistor between input node (input) and the 
    # inverting node (in_minus).
    circuit.R('2', 'in_minus', circuit.gnd, 1e3)    

    # Add an input sinewave voltage source for transient analysis.
    circuit.SinusoidalVoltageSource('V3', 
                                    'in_plus', circuit.gnd, 
                                    amplitude=1, 
                                    frequency=1e4, 
                                    offset=0)
    
    # Create a DC voltage source to simulate the negative supply.
    circuit.V('1', 'neg', circuit.gnd, 15)  # 15V DC source

    # Create a DC voltage source to simulate the positive supply.
    circuit.V('2', 'pos', circuit.gnd, -15)  # 15V DC source

    return circuit

def run_transient_analysis(circuit: Circuit) -> TransientAnalysis:
    try:
        # Initializes a simulator for the given circuit.
        simulator = circuit.simulator(temperature=25, nominal_temperature=25)

        # Set simulation parameters
        simulator.options(
            temp=25,      # sets the temperature for the simulation in degrees Celsius.
            trtol=1,      # sets the transient error tolerance.
            abstol=1e-12, # absolute tolerance for currents in the simulation.
            vntol=1e-6,   # specifies the voltage tolerance in the simulation.
            itl1=100,     # maximum number of iterations allowed for convergence 
                          # during DC analysis.
            itl4=50)      # maximum number of iterations allowed per timestep 
                          # during transient simulations.  

        # Perform an AC simulation of the circuit.
        simulation = simulator.transient(step_time=1e-6, end_time=1e-3)

        return simulation

    except Exception as e:
        print(f"Transient analysis error: {str(e)}")
        raise


def plot_bode_diagram(simulation):

    plt.figure(figsize=(8, 6))
    plt.title('Transient Analysis')

    # Extract the time and input voltage from analysis.
    time = np.array(simulation.time * 1e3)  # 1e3 is used to convert to ms.
    input_voltage = np.array(simulation['in_plus'])
    output_voltage = np.array(simulation['out'])

    # 
    plt.plot(time, input_voltage, label='Input Voltage')        
    plt.plot(time, output_voltage, label='Output Voltage')
    plt.grid(True, which="both", ls="-", alpha=0.6)
    plt.ylabel('Voltage (V)')
    plt.xlabel('Time (ms)')
    plt.tight_layout()
    plt.show()


def main():
    try:
        # Create the circuit
        circuit = create_circuit()

        # Run an AC analysis on the circuit
        simulation = run_transient_analysis(circuit)

        # Plot
        plot_bode_diagram(simulation)

    except Exception as e:
        print(f"Analysis error: {str(e)}")
        raise


if __name__ == "__main__":
    main()