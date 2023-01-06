#include <tests/common/SphereVerticesGenerator.h>

#include <glm/geometric.hpp>

#include <cmath>

namespace tests {
namespace common {
SphereVerticesGenerator::SphereVerticesGenerator(std::size_t slices, std::size_t stacks, float radius)
    : vertices(generateVertices(slices, stacks, radius))
    , normals(generateNormals(slices, stacks))
{
}
const double PI = acos(-1.0);

std::vector<glm::vec4> SphereVerticesGenerator::generateVertices(std::size_t slices, std::size_t stacks, float radius)
{
    const double sliceStep = PI * 2.0 / (double)slices;
    const double stackStep = PI / (double)stacks;
    std::vector<glm::vec4> result;

    auto pointOf = [radius](double theta, double phi) {
        return glm::vec4{radius * std::cos(theta) * std::cos(phi),
                         radius * std::cos(theta) * std::sin(phi),
                         radius * std::sin(theta), 1.0};
    };

    for (std::size_t stack = 0; stack < stacks; ++stack) {
        double theta = -PI / 2.0 + stack * stackStep;
        double thetaNext = -PI / 2.0 + (stack + 1) * stackStep;

        for (std::size_t slice = 0; slice < slices; ++slice) {
            double phi = slice * sliceStep;
            double phiNext = (slice + 1) * sliceStep;

            result.push_back(pointOf(theta, phi));
            result.push_back(pointOf(theta, phiNext));
            result.push_back(pointOf(thetaNext, phi));

            result.push_back(pointOf(thetaNext, phi));
            result.push_back(pointOf(theta, phiNext));
            result.push_back(pointOf(thetaNext, phiNext));
        }
    }

    return result;
}

std::vector<glm::vec4> SphereVerticesGenerator::generateNormals(std::size_t slices, std::size_t stacks)
{
    const double sliceStep = PI * 2.0 / (double)slices;
    const double stackStep = PI / (double)stacks;
    std::vector<glm::vec4> result;

    auto pointOf = [](double theta, double phi) {
        return glm::vec4{std::cos(theta) * std::cos(phi), //
                         std::cos(theta) * std::sin(phi), //
                         std::sin(theta), 1.0};
    };

    for (std::size_t stack = 0; stack < stacks; ++stack) {
        double theta = -PI / 2.0 + stack * stackStep;
        double thetaNext = -PI / 2.0 + (stack + 1) * stackStep;

        for (std::size_t slice = 0; slice < slices; ++slice) {
            double phi = slice * sliceStep;
            double phiNext = (slice + 1) * sliceStep;

            result.push_back(glm::normalize(pointOf(theta, phi)));
            result.push_back(glm::normalize(pointOf(theta, phiNext)));
            result.push_back(glm::normalize(pointOf(thetaNext, phi)));

            result.push_back(glm::normalize(pointOf(thetaNext, phi)));
            result.push_back(glm::normalize(pointOf(theta, phiNext)));
            result.push_back(glm::normalize(pointOf(thetaNext, phiNext)));
        }
    }

    return result;
}
}
}
