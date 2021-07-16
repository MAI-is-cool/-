// Основа универсального диффиренциального уравнения для последующего решения
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//*********************************************************
namespace StandartHelperLibrary.MathHelper
{
    /// <summary>
    /// Основа универсального диффиренциального уравнения для последующего решения 
    /// </summary>
    public abstract class TBaseDifferentialEquation
    {
        /// <summary>
        /// Крайняя левая точка диапазона "х" 
        /// </summary>
        public double Min_X { get; set; } = 0;
        /// <summary>
        /// Значение "у" в крайней левой точке 
        /// </summary>
        public double Min_Y { get; set; } = 1;
        /// <summary>
        /// Значение "z" в крайней левой точке 
        /// </summary>
        public double Min_Z { get; set; } = 1;
        /// <summary>
        /// Шаг сетки 
        /// </summary>
        public double Step { get; set; } = 0.1;
        /// <summary>
        /// Количество итераций c шагом Step  
        /// </summary>
        public int CountIterations { get; set; } = 10;           
        /// <summary>
        /// Округление до нужного знака, после запятой
        /// </summary>
        public int Rounding { get; set; } = 3;            
    }
}
