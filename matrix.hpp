#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cstddef>
#include <vector>
#include <ostream>
#include <iostream>

namespace dauphine
{
    class matrix
    {
    public:

        matrix(std::size_t nb_rows, std::size_t nb_cols);

    public:

        std::size_t nb_rows() const;
        std::size_t nb_cols() const;
        void resize(std::size_t nb_rows, std::size_t nb_cols);
        // Implementation: where is the const overload?
        double& operator()(std::size_t i, std::size_t j);
        // Design: should be a free function
        std::vector<double> produit_mat_vect(const std::vector<double>& v);
        

    private:

        std::size_t m_nb_rows;
        std::size_t m_nb_cols;
        std::vector<double> m_data;
    };

    

}


#endif
