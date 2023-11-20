

# LCMC (Long Chain Monte Carlo): Discover Alloy Ground State with Enhanced Monte-Carlo Simulation


## Table of Contents

- [Introduction](#Introduction)
- [Features](#Features)
- [Installation](#Installation)
- [Usage](#Usage)
- [Documentation](#Documentation)
- [Contributing](#Contributing)
- [License](#License)

## Introduction

Welcome to LCMC (Long Chain Monte Carlo), an innovative open-source package designed for finding alloy ground states using an enhanced Monte-Carlo simulation approach. Unlike traditional Monte Carlo simulations where each step is random, our implementation leverages the number of previous Monte Carlo steps to make more informed predictions, leading to significantly faster discovery of alloy ground states.

LCMC is a valuable tool for researchers and scientists in materials science and computational chemistry, providing an efficient and accurate way to explore and identify alloy configurations with minimal computational effort.

## Features

- Utilizes an enhanced Monte-Carlo simulation method to find alloy ground states.
- Incorporates previous Monte Carlo steps to make informed predictions for faster convergence.
- Benchmarked for improved efficiency compared to traditional Monte Carlo simulations.
- User-friendly interface and customizable parameters for tailored simulations.
- Open-source and well-documented for transparency and collaborative development.

## Installation

Regular installation

``` bash
make install
```

Developer installation
``` bash
make dev
```

## Usage

1. Create a "mcdft_parameters.json" file (you will find an example file on test folder).
2. 
``` bash
mcdft 
```


For comprehensive usage instructions, please consult the [User Guide](#).

## Documentation

Explore the comprehensive documentation to understand LCMC's features, usage, and customization options:

- [User Guide](#): Detailed instructions on how to effectively utilize LCMC.
- [API Reference](#): In-depth information about available functions and classes.
- [Examples](#): Real-world examples and use cases to help you get started.

## Contributing

We welcome contributions from the community! If you wish to contribute to the improvement of LCMC, please read our [Contributing Guidelines](#) to get started.

## License

LCMC is released under the [MIT License](#). Feel free to use, modify, and distribute it for your research or development projects.

---

Thank you for choosing LCMC for your materials science and computational chemistry research. We believe this package will accelerate your search for alloy ground states and provide valuable insights into your simulations. If you encounter any issues or have suggestions for improvement, please do not hesitate to [report them](#).

Happy simulating! 🧪🔬