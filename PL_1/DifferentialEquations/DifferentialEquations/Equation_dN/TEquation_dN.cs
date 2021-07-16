//Система дифференциальных уравнений с n переменными и n неизвестными
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//*********************************************************
namespace StandartHelperLibrary.MathHelper
{
    /// <summary>
    /// Система дифференциальных уравнений с n переменными и n неизвестными
    /// </summary>
    public class TEquation_dN : TBaseSystemDifferentialEquation, ISystemDifferentialEquation
    {
        /// <summary>
        /// Начальные значения Y
        /// </summary>
        public List<double> InitArray { get; set; } = new List<double>();
        /// <summary>
        /// Тело уравнения в виде делегата
        /// </summary>
        public AEquation_dN Equation { get; set; }
//-----------------------------------------------------------
        /// <summary>
        /// Рассчитать хранимое значение
        /// </summary>
        /// <param name="X">Переменная, по которой идет интегриривание</param>
        /// <param name="Y">Массив неизветных функций, зависящих от x </param>
        /// <returns>Рехзультат</returns>
        public double []  ComputeEquation(double X, double [] Y)
        {
            return Equation(X, Y);
        }

    }
}

