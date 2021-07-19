//
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//
using StandartHelperLibrary.MathHelper;
//***********************************************************
namespace PL_1
{
    class Program
    {
//----------------------------------------------------------
        static void Main(string[] args)
        {
            // Простой тест
            TDifferentialSolver.Debug(TDifferentialSolver.Example_dN_Residual());
            //TDifferentialSolverTwo.Debug(TDifferentialSolverTwo.Example_dXdYdZ());
            //
            Console.WriteLine("over");
            Console.ReadKey();
        }
//----------------------------------------------------------
    }
}
