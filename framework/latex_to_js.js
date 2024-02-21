#!/usr/bin/env node
const fs = require('fs');

// 3 because 0 is node path, 1 is script path and 2 is input file
if (process.argv.length != 3) {
    console.log(`ERROR: Please provide a latex file!\nUsage: ./latex_to_js example.tex`);
    process.exit(1);
}

console.log(`Reading file: ${process.argv[2]}`);
const latexString = fs.readFileSync(process.argv[2], 'utf8');

/*const latexString = `
\\begin{math}
    s(q)=\\frac{3}{2*\\pi}*
    \\begin{cases}
        \\frac{2}{3}-q^2+\\frac{1}{2}*q^3 & 0 \\leq q < 1\\\\
        \\frac{1}{6}*(2-q)^3 & 1 \\leq q < 2\\\\
        0 & q \\geq 2
    \\end{cases}
\\end{math}

\\begin{math}
    W(x1,x2,h)=\\frac{1}{h^3}*s(\\frac{||x1 - x2||_{2}}{h})
\\end{math}

\\begin{math}
    sgrad(q)=
    \\begin{cases}
            (\\frac{3}{2*\\pi}) * (-2*q+\\frac{3}{2}*q^2) & 0 \\leq q < 1\\\\
            (\\frac{3}{2*\\pi}) * (\\frac{-1}{2}*(2-q)^2) & 1 \\leq q < 2\\\\
            0 & q \\geq 2
    \\end{cases}
\\end{math}

\\begin{math}
    Wgrad(x1,x2,h)=
    \\begin{cases}
        \\begin{pmatrix} 0 \\\\ 0 \\\\ 0 \\end{pmatrix} & x1 = x2\\\\
        \\begin{pmatrix}
            \\frac{\\frac{1}{h^4}*sgrad(\\frac{||x1 - x2||_{2}}{h})}{||x1 - x2||_{2}} * (x1.x - x2.x) \\\\
            \\frac{\\frac{1}{h^4}*sgrad(\\frac{||x1 - x2||_{2}}{h})}{||x1 - x2||_{2}} * (x1.y - x2.y) \\\\
            \\frac{\\frac{1}{h^4}*sgrad(\\frac{||x1 - x2||_{2}}{h})}{||x1 - x2||_{2}} * (x1.z - x2.z)
        \\end{pmatrix}
        & x1 \\neq x2
    \\end{cases}
\\end{math}
`;*/

const mathBlockMatches = Array.from(latexString.matchAll(/\\begin{math}[\s\S\n]*?\\end{math}/gm), m => m[0]);

let resultingCode = '';

for (let mathBlock of mathBlockMatches) {
    resultingCode += generateFunctionFromMathBlock(mathBlock) + '\n';
}

// Custom functions (e.g. vector math) to make everything work!
resultingCode = getCustomCode() + resultingCode;

console.log(resultingCode);
fs.writeFileSync('output.cpp', resultingCode);

function generateFunctionFromMathBlock(mathBlock) {
    // Make sure there is exactly one function definition
    const functionDefinitions = Array.from(mathBlock.matchAll(/([A-Za-z0-9]*\s*\(.*\)\s*)=/gm), m => m[0]);
    if (functionDefinitions.length != 1) {
        console.log(`ERROR: Did not get exactly one function definition, be sure to include exactly one function in each math block!\nSupported format:\n\\begin{math}\nf(...)=...\n\\end{math}`);
        console.log(`\nMath Block:`);
        console.log(mathBlock);
        process.exit(1);
    }

    const bodyMatches = [...mathBlock.matchAll(/([A-Za-z0-9]*\s*\(.*\)\s*)=([\s\S]*)\\end{math}/gm)][0];
    if (bodyMatches[2].trim().replace(/\n*/g, '').length == 0) {
        console.log(`ERROR: Function misses definition!\nSupported format:\n\\begin{math}\nf(...)=...\n\\end{math}`);
        console.log(`\nMath Block:`);
        console.log(mathBlock);
        process.exit(1);
    }
    
    const latexFunctionDeclaration = bodyMatches[1];
    const latexFunctionDefinition = bodyMatches[2];

    let code = '';

    code += getFunctionBarebone(latexFunctionDeclaration);
    const body = getFunctionBody(latexFunctionDefinition);
    code = code.replace('$BODY', body);

    return code;
}

// NOTE: Assumes that vector/matrix entries are separated by \\ (Only supports single column vectors for now)
function replacePMatrices(latex) {
    while (true) {
        const pmatrixMatches = [...latex.matchAll(/\\begin{pmatrix}([\s\S]*?)\\end{pmatrix}/gm)];
        if (pmatrixMatches.length == 0) {
            return latex;
        }
        const match = pmatrixMatches[0][1];
        const splitted = match.split('\\\\');
        if (splitted.length != 3) {
            console.log(`ERROR: Currently only column vectors with exactly 3 entries are supported!`);
            console.log('\nVector:');
            console.log(latex);
            process.exit(1);
        }
       
        const code = `{${splitted[0]}, ${splitted[1]}, ${splitted[2]}}`;
        latex = latex.replace(pmatrixMatches[0][0], code);
    }
}

function replaceCondition(latex) {
    while(latex.indexOf('\\leq') >= 0) {
        latex = latex.replace('\\leq', '<=');
    }
    while(latex.indexOf('\\geq') >= 0) {
        latex = latex.replace('\\geq', '>=');
    }
    while(latex.indexOf(' = ') >= 0) {
        latex = latex.replace(' = ', ' == ');
    }

    const splitted = latex.split(' ');

    if (splitted.length == 5) {
        // e.g. x <= y < z
        return `${splitted[0]} ${splitted[1]} ${splitted[2]} && ${splitted[2]} ${splitted[3]} ${splitted[4]}`;
    } else if (splitted.length == 3) {
        // e.g. x >= 2 or x = 2
        return `${splitted[0]} ${splitted[1]} ${splitted[2]}`;
    } else if (splitted.length == 1 && splitted[0] == 'else') {
        return '';
    } else {
        console.log(`ERROR: Unsupported condition. Must be in one of the following formats:`);
        console.log('\tx <= y < z (example, also works with >= etc.)');
        console.log('\tx >= 2 (example, also works with < etc.)');
        console.log('\telse');
        console.log(`\nCondition:`);
        console.log(latex);
        process.exit(1);
    }
}

function replaceCases(latex) {
    while (true) {
        const casesMatches = [...latex.matchAll(/\\begin\{cases\}([\s\S]*?)\\end\{cases\}/gm)];
        if (casesMatches.length == 0) {
            return latex;
        }
        const caseMatch = casesMatches[0].map(m => m.trim())[1];
        const cases = caseMatch.split('\\\\\n');

        let caseCodes = [];

        for (const c of cases) {
            const splitted = c.split('&').map(s => s.trim());
            const formular = splitted[0];
            const condition = replaceCondition(splitted[1]);
            
            // NOTE: Last case is considered the 'else' case (without condition)
            if (c == cases[cases.length - 1]) {
                caseCodes.push(`${formular}`);
            } else {
                caseCodes.push(`((${condition}) ? ${formular}`);
            }
        }

        let code = caseCodes.join(' : ') + ')'.repeat(caseCodes.length - 1);
        latex = latex.replace(casesMatches[0][0], `${code}`);
    }
}

function replaceFractions(latex) {
    while (true) {
        const fractionMatches = [...latex.matchAll(/\\frac{(.*?)}{(.*?)}/gm)];
        if (fractionMatches.length == 0) {
            return latex;
        }
        const match = fractionMatches[0].map(m => m.trim());

        const topPart = match[1];
        const bottomPart = match[2];

        latex = latex.replace(match[0], `((${topPart})/(${bottomPart}))`);
    }
}

function replaceVectorNorm(latex) {
    while (true) {
        const vectorNormMatches = [...latex.matchAll(/\|\|(.*?)-(.*?)\|\|_\{2\}/gm)];
        if (vectorNormMatches.length == 0) {
            return latex;
        }
        const match = vectorNormMatches[0].map(m => m.trim());

        latex = latex.replace(match[0], `(distanceNorm(${match[1]},${match[2]}))`);
    }
}

function replacePI(latex) {
    while(latex.indexOf('\\pi') >= 0) {
        latex = latex.replace('\\pi', '(M_PI)');
    }
    return latex;
}

function replacePows(latex) {
    while(true) {
        const powMatches = [...latex.matchAll(/([A-Za-z0-9]+?)\^([0-9]+)/gm)];
        if (powMatches.length == 0) {
            return latex;
        }
        const match = powMatches[0];

        latex = latex.replace(match[0], `std::pow(${match[1]}, ${match[2]})`);
    }
}

function getFunctionBody(definition) {
    let code = definition;
  
    code = replacePMatrices(code);
    code = replaceCases(code);
    code = replaceVectorNorm(code);
    code = replaceFractions(code);
    code = replacePI(code);
    code = replacePows(code);

    if (code.endsWith('\n')) {
        code = code.substring(0, code.length - 1);
    }

    code = `    return ${code.trim()};`;

    // Remove newlines after return that would make the program incorrect
    code = code.replace(/return\s*\n/gm, 'return ');

    return code;
}

function getFunctionBarebone(declaration) {
    const openParantethesIndex = declaration.indexOf('(');
    const closeParantethesIndex = declaration.indexOf(')');

    const name = declaration.substring(0, openParantethesIndex);
    parameters = declaration.substring(openParantethesIndex + 1, closeParantethesIndex).split(',').map(param => param.trim());

    // NOTE: Hack to set correct parameter type
    parameters = parameters.map(param => { 
        if (param == 'x1' || param == 'x2' || param == 'v1' || param == 'v2' || param == 'n1' || param == 'n2') {
            return `Eigen::Vector3d ${param}`;
        } else {
            return `double ${param}`;
        }
    });
    console.log(parameters);

    return `double ${name}(${parameters.join(', ')}) {\n$BODY\n}`;
}

function getCustomCode() {
    return `\
#include <cmath>
#include <Eigen/Dense>

double distanceNorm(Eigen::Vector3d x1, Eigen::Vector3d x2) {
    return (x1-x2).norm();
}

`;
}
