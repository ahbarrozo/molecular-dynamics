#include "vector3.h"
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

int main() {
    std::cout << "=== Testing Vector3 ===\n\n";
    
    try {
        // Test 1: Create tensors with specific values
        std::cout << "Creating Vector3 in distinct ways...\n";
        
        Vector3 v1;
        Vector3 v2(1.0f, 2.0f, 3.0f);
        Vector3 v3(std::vector<float>{2.0f, 2.0f, 2.0f});
        
        std::cout << "Vector 1:\n";
        v1.print();
        
        std::cout << "Vector 2:\n";
        v2.print();
        
        std::cout << "Vector 3:\n";
        v3.print();
        
        //Test 2a: Test addition
        std::cout << "Testing addition...\n";
        Vector3 result = v1 + v2;
        
        std::cout << "Result (v1 + v2):\n";
        result.print();

        std::cout << "Verifying results...\n";
        
        assert(abs(result.x() - 1.0f) < 1e-6);
        assert(abs(result.y() - 2.0f) < 1e-6);  
        assert(abs(result.z() - 3.0f) < 1e-6);

        std::cout << "Testing in-place addition using previous results...\n";
        result += v3;
        
        std::cout << "Result ((v1 + v2) += v3):\n";
        result.print();

        std::cout << "Verifying results...\n";
        
        assert(abs(result.x() - 3.0f) < 1e-6);
        assert(abs(result.y() - 4.0f) < 1e-6);  
        assert(abs(result.z() - 5.0f) < 1e-6);

        std::cout << "All addition tests passed!\n\n";

        // Test 2b: Test subtraction
        std::cout << "Testing subtraction...\n";
        result = v1 - v2;
        
        std::cout << "Result (v1 - v2):\n";
        result.print();
        
        // Test 3: Verify results manually
        std::cout << "Verifying results...\n";
        
        assert(abs(result.x() + 1.0f) < 1e-6);
        assert(abs(result.y() + 2.0f) < 1e-6);  
        assert(abs(result.z() + 3.0f) < 1e-6);
        
        std::cout << "Testing in-place subtraction using previous results...\n";
        result -= v3;
        
        std::cout << "Result ((v1 - v2) -= v3):\n";
        result.print();

        std::cout << "Verifying results...\n";

        assert(abs(result.x() + 3.0f) < 1e-6);
        assert(abs(result.y() + 4.0f) < 1e-6);  
        assert(abs(result.z() + 5.0f) < 1e-6);

        std::cout << "All subtraction tests passed!\n\n";

        // Test 2c: Test scalar multiplication
        std::cout << "Testing scalar multiplication...\n";
        result = v2 * 2.5f;
        
        std::cout << "Result (v2 * 2.5):\n";
        result.print();
        
        // Test 3: Verify results manually
        std::cout << "Verifying results...\n";
        
        assert(abs(result.x() - 2.5f) < 1e-6);
        assert(abs(result.y() - 5.0f) < 1e-6);  
        assert(abs(result.z() - 7.5f) < 1e-6);
        
        std::cout << "Testing scalar multiplication commutation...\n";
        result = 2.5f * v2;
        
        std::cout << "Result (2.5 * t1 == t1 * 2.5):\n";

        assert(abs(result.x() - 2.5f) < 1e-6);
        assert(abs(result.y() - 5.0f) < 1e-6);  
        assert(abs(result.z() - 7.5f) < 1e-6);

        std::cout << "Testing in-place scalar multiplication using previous results..\n";
        result *= 2.0f;
        
        std::cout << "Result ((2.5 * t1) * 2):\n";
        result.print();

        assert(abs(result.x() - 5.0f) < 1e-6);
        assert(abs(result.y() - 10.0f) < 1e-6);  
        assert(abs(result.z() - 15.0f) < 1e-6);

        std::cout << "All scalar multiplication tests passed!\n\n";

        // Test 2d: Test vector products
        std::cout << "Testing dot products...\n";
        float dot_prod;
        dot_prod = v2.dot(v3);
        
        std::cout << "Result (v1 . v2): " << dot_prod << "\n";
        
        // Test 3: Verify results manually
        std::cout << "Verifying results...\n";
        
        assert(abs(dot_prod - 12.0f) < 1e-6);

        std::cout << "Testing cross product...\n";
        result = v2.cross(v3);
        
        std::cout << "Result (v2 x v3):\n";
        result.print();
        
        // Test 3: Verify results manually
        std::cout << "Verifying results...\n";

        assert(abs(result.x() + 2.0f) < 1e-6);
        assert(abs(result.y() - 4.0f) < 1e-6);  
        assert(abs(result.z() + 2.0f) < 1e-6); 
        
        std::cout << "All vector products tests passed!\n\n";

        // Test 4: Test with different shapes
        std::cout << "Testing PBC for cubic box of size 2...\n";

        Vector3 box(2.0f, 2.0f, 2.0f);

        result = v1.PBCWrap(box);
        std::cout << "Vector 1 after PBC: ";
        result.print();

        result = v2.PBCWrap(box);
        std::cout << "Vector 2 after PBC: ";
        result.print();
        result = v3.PBCWrap(box);
        std::cout << "Vector 3 after PBC: ";
        result.print();

        // std::cout << "Vector 2: ";
        // vec2.print();
        // std::cout << "Sum: ";
        // vec_result.print();
        
        // // Verify: [1+4, 2+5, 3+6] = [5, 7, 9]
        // assert(abs(vec_result({0}) - 5.0f) < 1e-6);
        // assert(abs(vec_result({1}) - 7.0f) < 1e-6);
        // assert(abs(vec_result({2}) - 9.0f) < 1e-6);
        
        // std::cout << "✓ Vector addition tests passed!\n\n";

        // // Test matrix multiplication
        // std::cout << "Test: Matrix multiplication...\n";
        // Tensor A({1.0f, 2.0f, 3.0f, 4.0f}, {2, 2});
        // Tensor B({5.0f, 6.0f, 7.0f, 8.0f}, {2, 2});

        // std::std::cout << "Matrix A:\n"; A.print();
        // std::std::cout << "Matrix B:\n"; B.print();

        // Tensor C = A.matrix_multiplication(B);
        // std::std::cout << "A * B:\n"; C.print();

        // // Verify results
        // assert(std::abs(C({0, 0}) - 19.0f) < 1e-6);  // 1*5 + 2*7 = 19
        // assert(std::abs(C({0, 1}) - 22.0f) < 1e-6);  // 1*6 + 2*8 = 22
        // assert(std::abs(C({1, 0}) - 43.0f) < 1e-6);  // 3*5 + 4*7 = 43
        // assert(std::abs(C({1, 1}) - 50.0f) < 1e-6);  // 3*6 + 4*8 = 50

        // std::std::cout << "✓ Matrix multiplication test passed!\n";
        
        std::cout << "All tests successful!\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}