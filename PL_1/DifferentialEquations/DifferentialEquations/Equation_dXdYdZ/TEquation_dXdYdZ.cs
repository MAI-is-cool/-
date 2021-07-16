// Дифференциальное уравнение вида d^2y/dx^2+p*dy/dx + q = f(x,y) 
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//*********************************************************
namespace StandartHelperLibrary.MathHelper
{
    /// <summary>
    /// Дифференциальное уравнение вида dd^2y/dx^2+p*dy/dx + q = f(x,y) 
    /// </summary>
    public class TEquation_dXdYdZ : TBaseDifferentialEquation, IDifferentialEquation
    {
        /// <summary>
        /// Тело уравнения в виде делегата
        /// </summary>
        public AEquation_dXdYdZ Equation { get; set; }
        //-----------------------------------------------------------
        /// <summary>
        /// Рассчитать хранимое уравнение
        /// </summary>
        /// <param name="Values">Значения</param>
        /// <returns>Результат</returns>
        public double ComputeEquation(params double[] Values)
        {
            return Equation(Values[0], Values[1], Values[2]);
        }

    }
}

