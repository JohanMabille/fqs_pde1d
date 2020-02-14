#include "matrix.hpp"
#include <ostream>

namespace dauphine
{
    matrix::matrix(std::size_t nb_rows, std::size_t nb_cols)
        : m_nb_rows(nb_rows),
          m_nb_cols(nb_cols),
          m_data(nb_rows * nb_cols)
    {
    }
    std::size_t matrix::nb_rows() const
    {
        return m_nb_rows;
    }

    std::size_t matrix::nb_cols() const
    {
        return m_nb_cols;
    }
    void matrix::resize(std::size_t nb_rows,std::size_t nb_cols)
    {
        m_nb_rows=nb_rows;
        m_nb_cols=nb_cols;
        
    }

    std::ostream& operator<<(std::ostream& out, matrix& m)
    {
        for(std::size_t i = 0; i < m.nb_rows(); ++i)
        {
            for(std::size_t j = 0; j < m.nb_cols(); ++j)
            {
                out << m(i, j) << ", ";
            }
            out << std::endl;
        }
        return out;
    }
    double& matrix::operator()(std::size_t i,std::size_t j)
    {
        return  m_data[i*m_nb_rows+j];
    }

    // Implementation: this has not been declared in the header
    // and therefore it cannot be used
    matrix operator+(double lhs, matrix& m)
    {
        for(std::size_t i = 0; i < m.nb_rows(); ++i)
        {
            for(std::size_t j = 0; j < m.nb_cols(); ++j)
            {
                m(i, j)=m(i,j)+lhs;
            }
        }
        return m;
    }

    std::vector<double> matrix::produit_mat_vect(const std::vector<double>& v)
    {
    if (v.size()==m_nb_cols)
        {
        std::vector<double> vect(m_nb_rows);
        for(std::size_t i=0;i<m_nb_rows;++i)
        {
            double sum=0;
            for(std::size_t j=0;j<m_nb_cols;++j)
            {
                sum+=m_data[i * m_nb_cols + j]*v[j];
            }
            vect[i]=sum;
        }
        return vect;
        }
    else
        {
            std::cout<< "Error of m_nb_cols dimensions";
        }
    }

}







