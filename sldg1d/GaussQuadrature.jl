# Funtion of Gauss Legendre quadrature
module GaussQuadrature

using Todo

## the points of Gauss Legendre quadrature nodes in [-1/2, 1/2]
todo"仿照FastGaussQuadrature 将我们的代码改写成更加灵活的的形式"
function gaussLegendreQuadrature(p::Integer)
    gaussLegendreMatrix = zeros(p + 1, p + 1)
    if p == 0
        gaussLegendreMatrix[1, 1] = 0.0
        gaussLegendreMatrix[1, 2] = 1.0

        return gaussLegendreMatrix
    elseif p == 1
        gaussLegendreMatrix[1, 1] = -sqrt(1.0/3.0) * 0.5
        gaussLegendreMatrix[2, 1] = sqrt(1.0/3.0) * 0.5

        gaussLegendreMatrix[1, 2] = 0.5
        gaussLegendreMatrix[2, 2] = 0.5

        return gaussLegendreMatrix
    elseif p == 2
        gaussLegendreMatrix[1, 1] = -sqrt(3.0/5.0) * 0.5
        gaussLegendreMatrix[2, 1] = 0.0
        gaussLegendreMatrix[3, 1] = sqrt(3.0/5.0) * 0.5

        gaussLegendreMatrix[1, 2] = 5.0/18.0
        gaussLegendreMatrix[2, 2] = 4.0/9.0        
        gaussLegendreMatrix[3, 2] = 5.0/18.0

        return gaussLegendreMatrix
    elseif p == 3
        gaussLegendreMatrix[1, 1] = -sqrt(3.0/7.0 + 2.0 / 7.0 * sqrt(1.2)) * 0.5
        gaussLegendreMatrix[2, 1] = -sqrt(3.0/7.0 - 2.0 / 7.0 * sqrt(1.2)) * 0.5
        gaussLegendreMatrix[3, 1] = sqrt(3.0/7.0 - 2.0 / 7.0 * sqrt(1.2)) * 0.5
        gaussLegendreMatrix[4, 1] = sqrt(3.0/7.0 + 2.0 / 7.0 * sqrt(1.2)) * 0.5

        gaussLegendreMatrix[1, 2] = (18.0 - sqrt(30.0)) / 36.0 * 0.5
        gaussLegendreMatrix[2, 2] = (18.0 + sqrt(30.0)) / 36.0 * 0.5
        gaussLegendreMatrix[3, 2] = (18.0 + sqrt(30.0)) / 36.0 * 0.5
        gaussLegendreMatrix[4, 2] = (18.0 - sqrt(30.0)) / 36.0 * 0.5

        return gaussLegendreMatrix
    elseif p == 4
        gaussLegendreMatrix[1, 1] = -1.0/3.0 * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0))*0.5
        gaussLegendreMatrix[2, 1] = -1.0/3.0 * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0))*0.5
        gaussLegendreMatrix[3, 1] = 0.0
        gaussLegendreMatrix[4, 1] = 1.0/3.0 * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0))*0.5
        gaussLegendreMatrix[5, 1] = 1.0/3.0 * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0))*0.5

        gaussLegendreMatrix[1, 2] = (322.0 - 13.0 * sqrt(70.0)) / 900.0 * 0.5
        gaussLegendreMatrix[2, 2] = (322.0 + 13.0 * sqrt(70.0)) / 900.0 * 0.5
        gaussLegendreMatrix[3, 2] = 128.0 / 225.0 * 0.5
        gaussLegendreMatrix[4, 2] = (322.0 + 13.0 * sqrt(70.0)) / 900.0 * 0.5
        gaussLegendreMatrix[5, 2] = (322.0 - 13.0 * sqrt(70.0)) / 900.0 * 0.5

        return gaussLegendreMatrix
    end
end

## the points of Gauss Lobatto quadrature nodes in [-1/2, 1/2]
todo"仿照FastGaussQuadrature 将我们的代码改写成更加灵活的的形式"
function gaussLobattoQuadrature(p::Integer)
    gaussLobattoMatrix = zeros(p + 1, p + 1)
    if p == 1
        gaussLobattoMatrix[1,1] = -0.5
        gaussLobattoMatrix[2,1] = 0.5

        gaussLobattoMatrix[1,2] = 0.5
        gaussLobattoMatrix[2,2] = 0.5

        return gaussLobattoMatrix
    elseif p == 2
        gaussLobattoMatrix[1,1] = -0.5
        gaussLobattoMatrix[2,1] = 0.0
        gaussLobattoMatrix[3,1] = 0.5

        gaussLobattoMatrix[1,2] = 1.0/6.0
        gaussLobattoMatrix[2,2] = 2.0/3.0
        gaussLobattoMatrix[3,2] = 1.0/6.0
        return gaussLobattoMatrix
    elseif p == 3
        gaussLobattoMatrix[1,1] = -0.5
        gaussLobattoMatrix[2,1] = -sqrt(5.0)/10.0
        gaussLobattoMatrix[3,1] = sqrt(5.0)/10.0
        gaussLobattoMatrix[4,1] = 0.5

        gaussLobattoMatrix[1,2] = 1.0/12.0
        gaussLobattoMatrix[2,2] = 5.0/12.0
        gaussLobattoMatrix[3,2] = 5.0/12.0
        gaussLobattoMatrix[4,2] = 1.0/12.0
        return gaussLobattoMatrix
    elseif p == 4
        gaussLobattoMatrix[1,1] = -0.5
        gaussLobattoMatrix[2,1] = -sqrt(21.0)/14.0
        gaussLobattoMatrix[3,1] = 0.0
        gaussLobattoMatrix[4,1] = sqrt(21.0)/14.0
        gaussLobattoMatrix[5,1] = 0.5

        gaussLobattoMatrix[1,2] = 1.0/20.0
        gaussLobattoMatrix[2,2] = 49.0/180.0
        gaussLobattoMatrix[3,2] = 64.0/180.0
        gaussLobattoMatrix[4,2] = 49.0/180.0
        gaussLobattoMatrix[5,2] = 1.0/20.0
        return gaussLobattoMatrix
    end
end

end